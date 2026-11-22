// Solid-wall adapter for the v2112 libAcoustics Farassat 1A implementation.
#include "FfowcsWilliamsHawkings.H"
#include "addToRunTimeSelectionTable.H"
#include "wallPolyPatch.H"
#include "quaternion.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include <cmath>
#include <limits>

namespace Foam
{
namespace functionObjects
{

class hbSolidWallFWH : public FfowcsWilliamsHawkings
{
    dictionary settings_;

    static dictionary backendSettings(const dictionary& settings)
    {
        const wordReList patches(settings.lookup("wallPatches"));
        if (patches.empty() || settings.subDict("observers").empty())
            FatalErrorInFunction << "Specify wallPatches and observers" << exit(FatalError);
        if (mag(settings.getOrDefault<vector>("U0", Zero)) > SMALL)
            FatalErrorInFunction << "This solid-wall Farassat 1A replay requires U0 = 0"
                << exit(FatalError);
        for (const word key : {word("rhoInf"), word("c0")})
        {
            const scalar value = settings.get<scalar>(key);
            if (!std::isfinite(value) || value <= 0)
                FatalErrorInFunction << key << " must be positive" << exit(FatalError);
        }
        // The upstream constructor sizes its histories from a sampledPatch.
        // The replay supplies every face in a single consistent material order.
        OStringStream text;
        text << "patches (); surfaces (walls { type patch; patches " << patches
             << "; interpolate false; triangulate false; });";
        IStringStream stream(text.str());
        dictionary d(stream);
        d.add("region", word("subTimeLevel0"));
        d.add("rho", word("rhoInf"));
        d.add("rhoInf", settings.get<scalar>("rhoInf"));
        d.add("c0", settings.get<scalar>("c0"));
        d.add("CofR", vector::zero);
        d.add("pInf", settings.get<scalar>("pInf"));
        d.add("formulationType", word("Farassat1AFormulation"));
        d.add("interpolationScheme", word("cell"));
        d.add("fixedResponseDelay", false);
        d.add("timeStart", 0);
        d.add("timeEnd", GREAT);
        d.add("observers", settings.subDict("observers"));
        return d;
    }

    static tensor rotation(const vector& omega, scalar time)
    {
        const scalar speed = mag(omega);
        return speed > 0 ? quaternion(omega/speed, speed*time).R() : tensor::I;
    }

public:
    TypeName("hbSolidWallFWH");

    hbSolidWallFWH(const word& name, const Time& clock, const dictionary& dict)
    :
        FfowcsWilliamsHawkings(name, clock, backendSettings(dict)),
        settings_(dict)
    {}

    bool execute() override
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        Time& clock = const_cast<Time&>(mesh.time());
        const scalarList phases(settings_.lookup("phaseTimes"));
        const scalarList frequencies(settings_.lookup("angularFrequencies"));
        const List<scalarList> re(settings_.lookup("transformReal"));
        const List<scalarList> im(settings_.lookup("transformImag"));
        if (phases.size() < 2 || frequencies.empty()
         || re.size() != frequencies.size() || im.size() != frequencies.size())
            FatalErrorInFunction << "Invalid HB temporal transform" << exit(FatalError);
        forAll(frequencies, h)
            if (re[h].size() != phases.size() || im[h].size() != phases.size())
                FatalErrorInFunction << "Invalid HB transform row " << h << exit(FatalError);
        const scalar period = settings_.get<scalar>("period");
        const scalar offset = settings_.get<scalar>("geometryOffset");
        const scalar geomTime = phases.first() + offset;
        const vector origin(settings_.lookup("origin"));
        const vector omega(settings_.lookup("omega"));
        const label samples = settings_.get<label>("samplesPerPeriod");
        const label periods = settings_.get<label>("periods");
        const scalar tol = settings_.getOrDefault<scalar>("geometryTolerance", 1e-7);
        if (!std::isfinite(period) || period <= 0 || samples < 3 || periods < 1
         || tol <= 0 || !std::isfinite(tol) || !std::isfinite(mag(origin)+mag(omega)+pInf_))
            FatalErrorInFunction << "Invalid replay interval, geometry or pressure reference" << exit(FatalError);
        const scalar dt = period/samples;
        for (scalar w : frequencies)
        {
            const scalar cycles = w*period/constant::mathematical::twoPi;
            if (!std::isfinite(cycles) || mag(cycles-std::round(cycles)) > 1e-7
             || 2*mag(cycles) >= samples)
                FatalErrorInFunction << "Unresolved or nonperiodic HB frequency " << w << exit(FatalError);
        }
        const wordReList patterns(settings_.lookup("wallPatches"));
        const labelList ids(mesh.boundaryMesh().patchSet(patterns).sortedToc());
        label faces = 0;
        for (label id : ids)
        {
            const polyPatch& pp = mesh.boundaryMesh()[id];
            if (!isA<wallPolyPatch>(pp))
                FatalErrorInFunction << pp.name() << " is not a solid wall" << exit(FatalError);
            faces += pp.size();
        }
        if (!returnReduce(faces, sumOp<label>()) || faces != controlSurfaces_[0].Cf().size())
            FatalErrorInFunction << "Empty wall selection or inconsistent sampled face count" << exit(FatalError);
        vectorField radius(faces), area0(faces);
        List<scalarField> pressure(phases.size());
        forAll(phases, k)
        {
            const fvMesh& phase = clock.lookupObject<fvMesh>("subTimeLevel"+Foam::name(k));
            const volScalarField& p = phase.lookupObject<volScalarField>("p");
            if (p.dimensions() != dimPressure && p.dimensions() != sqr(dimVelocity))
                FatalErrorInFunction << "Unsupported p dimensions " << p.dimensions() << exit(FatalError);
            const scalar scale = p.dimensions() == dimPressure ? 1 : rhoRef_;
            pressure[k].setSize(faces);
            const tensor turn = rotation(omega, phases[k]-phases[0]);
            label f = 0;
            for (label id : ids)
            {
                const polyPatch& ref = mesh.boundaryMesh()[id];
                const label pid = phase.boundaryMesh().findPatchID(ref.name());
                if (pid < 0 || phase.boundaryMesh()[pid].size() != ref.size())
                    FatalErrorInFunction << "Patch mismatch " << ref.name() << exit(FatalError);
                const polyPatch& pp = phase.boundaryMesh()[pid];
                forAll(pp, j)
                {
                    if (!k)
                    {
                        radius[f] = pp.faceCentres()[j]-origin;
                        area0[f] = -pp.faceAreas()[j]; // solid into fluid
                    }
                    if (mag(area0[f]) <= VSMALL
                     || mag(pp.faceCentres()[j]-origin-(turn & radius[f])) > tol
                     || mag(pp.faceAreas()[j]+(turn & area0[f])) > 1e-6*mag(area0[f]))
                        FatalErrorInFunction << "Rigid-wall geometry mismatch at phase " << k
                            << ", patch " << ref.name() << ", face " << j << exit(FatalError);
                    pressure[k][f] = scale*p.boundaryField()[pid][j]-pInf_;
                    if (!std::isfinite(pressure[k][f]))
                        FatalErrorInFunction << "Non-finite wall pressure" << exit(FatalError);
                    ++f;
                }
            }
        }
        // Cache the retained pressure modes, including the steady material load.
        List<scalarField> a(frequencies.size()), b(frequencies.size());
        forAll(frequencies, h)
        {
            a[h].setSize(faces, 0); b[h].setSize(faces, 0);
            forAll(phases, k)
            {
                a[h] += re[h][k]*pressure[k];
                b[h] += im[h][k]*pressure[k];
            }
        }
        pressure.clear();
        scalar maxRadius = 0, maxSpeed = 0, minDelayDistance = GREAT;
        for (const vector& r : radius)
        {
            maxRadius = max(maxRadius, mag(r));
            maxSpeed = max(maxSpeed, mag(omega ^ r));
            for (const SoundObserver& obs : observers_)
            {
                const vector d = obs.position()-origin;
                scalar dmin = mag(d-r);
                if (mag(omega) > 0)
                {
                    const vector axis = omega/mag(omega);
                    dmin = Foam::sqrt(sqr((d-r)&axis) + sqr(mag(d^axis)-mag(r^axis)));
                }
                minDelayDistance = min(minDelayDistance, dmin);
            }
        }
        reduce(maxRadius, maxOp<scalar>()); reduce(maxSpeed, maxOp<scalar>());
        reduce(minDelayDistance, minOp<scalar>());
        if (maxSpeed >= c0_ || minDelayDistance <= SMALL)
            FatalErrorInFunction << "Requires subsonic walls and off-surface observers" << exit(FatalError);
        fwhFormulation& fw = fwhFormulationPtr_();
        scalar maxDelay = 0;
        forAll(observers_, i)
        {
            if (!std::isfinite(mag(observers_[i].position())))
                FatalErrorInFunction << "Non-finite observer coordinates" << exit(FatalError);
            fw.tauMax_[i] = (mag(observers_[i].position()-origin)+maxRadius)/c0_;
            maxDelay = max(maxDelay, fw.tauMax_[i]);
        }
        // Discard the propagation transient and two source-derivative startup steps.
        const scalar warmPeriods = std::ceil((maxDelay+3*dt)/period);
        if (!std::isfinite(warmPeriods) || (warmPeriods+periods)*samples >= std::numeric_limits<label>::max())
            FatalErrorInFunction << "Acoustic sample count exceeds label range" << exit(FatalError);
        const label warm = label(warmPeriods)*samples;
        const label recordSize = periods*samples;
        label orders = 0;
        scalar bpf = 0, pRef = 0;
        if (const dictionary* tones = settings_.findDict("tones"))
        {
            const label blades = tones->get<label>("bladeCount");
            orders = tones->get<label>("orders");
            pRef = tones->get<scalar>("referencePressure");
            bpf = blades*mag(omega)/constant::mathematical::twoPi;
            const scalar cycles = bpf*dt*recordSize;
            if (blades < 1 || orders < 1 || pRef <= 0 || !std::isfinite(pRef)
             || bpf <= 0 || 2.0*orders*bpf*dt >= 1 || cycles < 1-SMALL
             || mag(cycles-std::round(cycles)) > 1e-8)
                FatalErrorInFunction << "Tones require resolved BPF orders, positive reference"
                    << " pressure and a whole number of blade passages" << exit(FatalError);
        }
        List<scalarField> cosine(observers_.size(), scalarField(orders, Zero));
        List<scalarField> sine(cosine);
        const scalar savedTime = clock.value();
        const label savedIndex = clock.timeIndex();
        const dimensionedScalar savedDt = clock.deltaT();
        clock.setDeltaT(dt, false);
        vectorField centres(faces), sf(faces);
        scalarField p(faces), density(faces, rhoRef_);
        vectorField& speed = vS_[0];
        if (Pstream::master()) analogyOutPtr_().precision(14);
        for (label step=0; step < warm + periods*samples; ++step)
        {
            const scalar t = step*dt;
            clock.setTime(t, step);
            const tensor turn = rotation(omega, t-geomTime);
            p = Zero;
            forAll(frequencies, h)
                p += Foam::cos(frequencies[h]*t)*a[h] - Foam::sin(frequencies[h]*t)*b[h];
            forAll(centres, f)
            {
                const vector r = turn & radius[f];
                centres[f] = origin+r;
                sf[f] = turn & area0[f];
                fw.ni_[0].value(f) = sf[f]/mag(sf[f]);
                speed[f] = omega ^ r;
            }
            if (step >= warm && Pstream::master()) analogyOutPtr_() << t;
            forAll(observers_, i)
            {
                forAll(centres, f)
                {
                    const vector r = observers_[i].position()-centres[f];
                    fw.robs_[i][0][f] = r;
                    fw.magrobs_[i][0][f] = mag(r);
                    fw.tobs_[i][0][f] = t+mag(r)/c0_;
                }
                const scalar value = fw.observerAcousticPressure(sf, speed, density, p, i, 0, t);
                if (!std::isfinite(value))
                    FatalErrorInFunction << "Non-finite observer pressure at " << t << exit(FatalError);
                if (step >= warm && Pstream::master())
                {
                    analogyOutPtr_() << ' ' << value;
                    const scalar angle = constant::mathematical::twoPi*bpf*(step-warm)*dt;
                    for (label h=0; h<orders; ++h)
                    {
                        cosine[i][h] += value*Foam::cos((h+1)*angle);
                        sine[i][h] -= value*Foam::sin((h+1)*angle);
                    }
                }
            }
            if (step >= warm && Pstream::master()) analogyOutPtr_() << nl;
            fw.clearExpiredData();
        }
        clock.setTime(savedTime, savedIndex);
        clock.setDeltaT(savedDt, false);
        if (Pstream::master()) analogyOutPtr_().flush();
        if (orders && Pstream::master())
        {
            OFstream out(analogyOutPtr_().name().path()/(name()+"-bpf.dat"));
            out.precision(14);
            out << "observer order frequency_Hz rms_Pa SPL_dB phase_deg_from_record_start reference_Pa" << nl;
            forAll(observers_, i)
                for (label h=0; h<orders; ++h)
                {
                    const scalar rms = std::hypot(cosine[i][h], sine[i][h])*std::sqrt(2.)/recordSize;
                    const scalar db = rms > 0 ? 20*std::log10(rms/pRef)
                        : -std::numeric_limits<scalar>::infinity();
                    const scalar phase = rms > 0 ? std::atan2(sine[i][h], cosine[i][h])
                        *180/constant::mathematical::pi : std::numeric_limits<scalar>::quiet_NaN();
                    out << observers_[i].name() << ' ' << h+1 << ' ' << (h+1)*bpf << ' '
                        << rms << ' ' << db << ' ' << phase << ' ' << pRef << nl;
                }
        }
        Info<< "HB acoustics wrote " << periods*samples << " observer samples after "
            << warm << " warm-up samples. dt = " << dt << " s" << endl;
        return true;
    }
};

defineTypeNameAndDebug(hbSolidWallFWH, 0);
addToRunTimeSelectionTable(functionObject, hbSolidWallFWH, dictionary);

} // namespace functionObjects
} // namespace Foam
