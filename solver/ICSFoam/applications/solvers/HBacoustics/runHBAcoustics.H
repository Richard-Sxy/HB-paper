#ifndef runHBAcoustics_H
#define runHBAcoustics_H

#include "functionObject.H"
#include "IOHBZoneList.H"
#include "IOdictionary.H"

namespace Foam
{

inline void runHBAcoustics
(
    Time& clock, const IOHBZoneList& hb, scalar motionOffset, bool required = false
)
{
    IOobject io("hbAcousticsDict", clock.system(), clock, IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE, false);
    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        if (required)
            FatalErrorInFunction << "Replay requires system/hbAcousticsDict" << exit(FatalError);
        return;
    }
    IOdictionary settings(io);
    if (!settings.getOrDefault<Switch>("enabled", true))
    {
        if (required)
            FatalErrorInFunction << "Replay requires enabled true" << exit(FatalError);
        return;
    }

    const HBZoneList& zones = hb;
    const HBZone* selected = nullptr;
    forAll(zones, i)
    {
        if (!zones[i].active()) continue;
        if (selected || !hb.subDict(zones[i].name()).getOrDefault<bool>("allMesh", false))
            FatalErrorInFunction << "Acoustic replay requires one active whole-mesh HB zone"
                << exit(FatalError);
        selected = &zones[i];
    }
    if (!selected) FatalErrorInFunction << "No active HB zone" << exit(FatalError);
    const auto& phase = hb.selectedSnapshots();
    const auto& e = selected->E();
    const auto& freq = selected->omegaList();
    if (phase.size() < 2 || e.m() != freq.size() || e.n() != phase.size())
        FatalErrorInFunction << "Inconsistent HB transform dimensions" << exit(FatalError);

    dictionary config(settings);
    config.set("type", word("hbSolidWallFWH"));
    config.set("libs", wordList(1, "libAcoustics.so"));
    config.set("region", word("subTimeLevel0"));
    config.set("phaseTimes", phase);
    config.set("angularFrequencies", freq);
    config.set("period", phase.size()*(phase[1]-phase[0]));
    config.set("geometryOffset", motionOffset);
    List<scalarList> re(freq.size()), im(freq.size());
    forAll(freq, h)
    {
        re[h].setSize(phase.size());
        im[h].setSize(phase.size());
        forAll(phase, k)
        {
            re[h][k] = e[h][k].Re();
            im[h][k] = e[h][k].Im();
        }
    }
    config.set("transformReal", re);
    config.set("transformImag", im);
    // functionObject::New does not load extra libraries if the type is already registered.
    const word name = config.getOrDefault<word>("outputName", "hbFwh_"+clock.timeName());
    autoPtr<functionObject> acoustic = functionObject::New(name, clock, config);
    if (!acoustic->execute())
        FatalErrorInFunction << "HB acoustic replay failed" << exit(FatalError);
}

}
#endif
