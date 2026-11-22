// Standalone HB and Jacobian-free Newton/GMRES algebra. C++11, no dependencies.
// All components supplied by the caller must already be nondimensionalized.
#ifndef IHB_HB_ALGEBRA_H
#define IHB_HB_ALGEBRA_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ihb
{
using Vec = std::vector<double>;
using Residual = std::function<Vec(const Vec&)>;
using Reporter = std::function<void(int, double, double)>;

inline bool finite(const Vec& a)
{
    for (double v : a) if (!std::isfinite(v)) return false;
    return true;
}

inline double dot(const Vec& a, const Vec& b)
{
    if (a.size() != b.size()) throw std::invalid_argument("dot: size mismatch");
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i]*b[i];
    return s;
}

inline double scaledNorm(const Vec& a, bool mean)
{
    // BLAS-style scaling avoids squaring overflow/underflow.
    double scale = 0.0, sum = 1.0;
    for (double value : a)
    {
        const double av = std::abs(value);
        if (!std::isfinite(av)) return av;
        if (av == 0.0) continue;
        if (av > scale)
        {
            const double ratio = scale/av;
            sum = 1.0+sum*ratio*ratio;
            scale = av;
        }
        else
        {
            const double ratio = av/scale;
            sum += ratio*ratio;
        }
    }
    if (scale == 0.0) return 0.0;
    return scale*std::sqrt(mean ? sum/double(a.size()) : sum);
}
inline double norm(const Vec& a) { return scaledNorm(a, false); }
inline double rms(const Vec& a) { return scaledNorm(a, true); }

inline std::vector<Vec> differentiation(int harmonics, double omega)
{
    if (harmonics < 0 || harmonics > (std::numeric_limits<int>::max()-1)/2
        || !std::isfinite(omega) || omega < 0.0)
        throw std::invalid_argument("HB requires nonnegative H and omega");
    const int n = 2*harmonics + 1;
    const double pi = std::acos(-1.0);
    std::vector<Vec> d(n, Vec(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = i+1; j < n; ++j)
        {
            double value = 0.0;
            for (int k = 1; k <= harmonics; ++k)
                value -= 2.0*omega*k*std::sin(2.0*pi*k*(i-j)/n)/n;
            d[i][j] = value;
            d[j][i] = -value;
        }
    return d;
}

struct Settings
{
    double absTol = 1e-9;
    double relTol = 1e-8;
    int maxNewton = 30;
    int gmresRestart = 60;
    int gmresMaxIterations = 300;
    double linearRelTol = 1e-3;
    double linearAbsTol = 1e-12;
    double fdStep = 1e-6;
    int maxLineSearch = 14;
    double armijo = 1e-4;
    double minStep = 1e-8;
};

struct Result
{
    bool converged = false;
    int iterations = 0;
    int linearIterations = 0;
    int evaluations = 0;
    double residual = std::numeric_limits<double>::infinity();
    std::string reason;
};

struct LinearResult
{
    Vec x;
    bool converged = false;
    int iterations = 0;
    double residual = std::numeric_limits<double>::infinity();
    std::string reason;
};

// Restarted, unpreconditioned GMRES. The state/residual scaling supplied by
// the caller acts as diagonal scaling; a physics block preconditioner can be
// added outside this routine. Residual values and tolerances are RMS norms.
inline LinearResult gmres
(
    const Residual& apply,
    const Vec& b,
    int restart,
    int maxIterations,
    double relTol,
    double absTol
)
{
    LinearResult out;
    const std::size_t n = b.size();
    out.x.assign(n, 0.0);
    if (n == 0 || !finite(b) || restart <= 0 || maxIterations <= 0
        || !std::isfinite(relTol) || !std::isfinite(absTol)
        || relTol < 0.0 || absTol < 0.0)
    {
        out.reason = "invalid GMRES input";
        return out;
    }
    const double sqrtN = std::sqrt(double(n));
    const double target = std::max(absTol, relTol*rms(b));
    Vec r = b;
    out.residual = rms(r);
    if (out.residual <= target)
    {
        out.converged = true;
        out.reason = "linear residual tolerance";
        return out;
    }

    while (out.iterations < maxIterations)
    {
        const int m = std::min
        (
            std::min(restart, int(n)), maxIterations-out.iterations
        );
        const double beta = norm(r);
        if (!std::isfinite(beta))
        {
            out.reason = "overflow in GMRES residual norm";
            return out;
        }
        std::vector<Vec> v(m+1, Vec(n, 0.0));
        std::vector<Vec> h(m+1, Vec(m, 0.0));
        Vec cs(m, 0.0), sn(m, 0.0), g(m+1, 0.0);
        for (std::size_t i = 0; i < n; ++i) v[0][i] = r[i]/beta;
        g[0] = beta;
        int used = 0;
        bool breakdown = false;
        for (int j = 0; j < m; ++j)
        {
            Vec w = apply(v[j]);
            ++out.iterations;
            if (w.size() != n || !finite(w))
            {
                out.reason = "nonfinite or mismatched GMRES operator";
                return out;
            }
            const double before = norm(w);
            if (!std::isfinite(before))
            {
                out.reason = "overflow in GMRES operator norm";
                return out;
            }
            // Two modified Gram-Schmidt passes limit loss of orthogonality.
            for (int pass = 0; pass < 2; ++pass)
                for (int k = 0; k <= j; ++k)
                {
                    const double projection = dot(w, v[k]);
                    h[k][j] += projection;
                    for (std::size_t i = 0; i < n; ++i)
                        w[i] -= projection*v[k][i];
                }
            h[j+1][j] = norm(w);
            breakdown = h[j+1][j] <=
                64.0*std::numeric_limits<double>::epsilon()*before;
            if (!breakdown && h[j+1][j] > 0.0)
                for (std::size_t i = 0; i < n; ++i)
                    v[j+1][i] = w[i]/h[j+1][j];
            else breakdown = true;
            for (int k = 0; k < j; ++k)
            {
                const double upper = cs[k]*h[k][j] + sn[k]*h[k+1][j];
                h[k+1][j] = -sn[k]*h[k][j] + cs[k]*h[k+1][j];
                h[k][j] = upper;
            }
            const double denominator = std::hypot(h[j][j], h[j+1][j]);
            if (!(denominator > 0.0) || !std::isfinite(denominator))
            {
                out.reason = "singular GMRES Arnoldi factor";
                return out;
            }
            cs[j] = h[j][j]/denominator;
            sn[j] = h[j+1][j]/denominator;
            h[j][j] = denominator;
            h[j+1][j] = 0.0;
            g[j+1] = -sn[j]*g[j];
            g[j] *= cs[j];
            used = j+1;
            if (std::abs(g[j+1])/sqrtN <= target || breakdown) break;
        }

        Vec y(used, 0.0);
        for (int j = used-1; j >= 0; --j)
        {
            double rhs = g[j];
            for (int k = j+1; k < used; ++k) rhs -= h[j][k]*y[k];
            if (h[j][j] == 0.0 || !std::isfinite(h[j][j]))
            {
                out.reason = "singular GMRES triangular factor";
                return out;
            }
            y[j] = rhs/h[j][j];
        }
        for (int j = 0; j < used; ++j)
            for (std::size_t i = 0; i < n; ++i)
                out.x[i] += y[j]*v[j][i];
        if (!finite(out.x))
        {
            out.reason = "nonfinite GMRES iterate";
            return out;
        }
        const Vec ax = apply(out.x);
        if (ax.size() != n || !finite(ax))
        {
            out.reason = "nonfinite or mismatched true linear residual";
            return out;
        }
        for (std::size_t i = 0; i < n; ++i) r[i] = b[i]-ax[i];
        out.residual = rms(r);
        // The Arnoldi estimate alone never certifies convergence.
        if (out.residual <= target)
        {
            out.converged = true;
            out.reason = "true linear residual tolerance";
            return out;
        }
        if (breakdown)
        {
            out.reason = "GMRES breakdown without true residual convergence";
            return out;
        }
    }
    out.reason = "GMRES iteration limit";
    return out;
}

inline bool valid(const Settings& s)
{
    return std::isfinite(s.absTol) && s.absTol >= 0.0
        && std::isfinite(s.relTol) && s.relTol >= 0.0
        && s.maxNewton >= 0 && s.gmresRestart > 0
        && s.gmresMaxIterations > 0
        && std::isfinite(s.linearRelTol) && s.linearRelTol > 0.0
        && s.linearRelTol < 1.0
        && std::isfinite(s.linearAbsTol) && s.linearAbsTol >= 0.0
        && std::isfinite(s.fdStep) && s.fdStep > 0.0
        && s.maxLineSearch > 0
        && std::isfinite(s.armijo) && s.armijo > 0.0 && s.armijo < 1.0
        && std::isfinite(s.minStep) && s.minStep > 0.0 && s.minStep <= 1.0;
}

inline Result newton
(
    Vec& x,
    const Residual& residual,
    const Settings& settings = Settings(),
    const Reporter& report = Reporter()
)
{
    Result out;
    if (x.empty() || !finite(x) || !valid(settings))
    {
        out.reason = "invalid Newton input";
        return out;
    }
    const auto evaluate = [&](const Vec& point) -> Vec
    {
        ++out.evaluations;
        return residual(point);
    };
    Vec f = evaluate(x);
    if (f.size() != x.size() || !finite(f))
    {
        out.reason = "nonfinite or mismatched initial residual";
        return out;
    }
    out.residual = rms(f);
    const double target = std::max
    (
        settings.absTol, settings.relTol*out.residual
    );
    if (report) report(0, out.residual, 0.0);
    for (int iteration = 0; iteration <= settings.maxNewton; ++iteration)
    {
        if (out.residual <= target)
        {
            out.converged = true;
            out.reason = "nonlinear residual tolerance";
            return out;
        }
        if (iteration == settings.maxNewton) break;

        const auto jv = [&](const Vec& v) -> Vec
        {
            const double vn = rms(v);
            if (vn == 0.0) return Vec(x.size(), 0.0);
            const double step = settings.fdStep*(1.0+rms(x))/vn;
            Vec xp = x, xm = x;
            for (std::size_t i = 0; i < x.size(); ++i)
            {
                xp[i] += step*v[i];
                xm[i] -= step*v[i];
            }
            const Vec fp = evaluate(xp), fm = evaluate(xm);
            if (fp.size() != x.size() || fm.size() != x.size()
                || !finite(fp) || !finite(fm))
                return Vec(x.size(), std::numeric_limits<double>::quiet_NaN());
            Vec result(x.size());
            for (std::size_t i = 0; i < x.size(); ++i)
                result[i] = (fp[i]-fm[i])/(2.0*step);
            return result;
        };
        Vec b = f;
        for (double& value : b) value = -value;
        const LinearResult linear = gmres
        (
            jv, b, settings.gmresRestart, settings.gmresMaxIterations,
            settings.linearRelTol, settings.linearAbsTol
        );
        out.linearIterations += linear.iterations;
        if (!linear.converged)
        {
            out.reason = "Newton linear solve failed: " + linear.reason;
            return out;
        }

        bool accepted = false;
        double alpha = 1.0;
        for (int trial = 0; trial < settings.maxLineSearch; ++trial)
        {
            if (alpha < settings.minStep) break;
            Vec candidate = x;
            for (std::size_t i = 0; i < x.size(); ++i)
                candidate[i] += alpha*linear.x[i];
            if (finite(candidate))
            {
                Vec trialF = evaluate(candidate);
                if (trialF.size() != x.size())
                {
                    out.reason = "mismatched line search residual";
                    return out;
                }
                const double trialNorm = finite(trialF)
                    ? rms(trialF) : std::numeric_limits<double>::infinity();
                if (trialNorm <= (1.0-settings.armijo*alpha)*out.residual
                    || trialNorm <= target)
                {
                    x.swap(candidate);
                    f.swap(trialF);
                    out.residual = trialNorm;
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if (!accepted)
        {
            out.reason = "Newton line search failed";
            return out;
        }
        out.iterations = iteration+1;
        if (report) report(out.iterations, out.residual, alpha);
    }
    out.reason = "Newton iteration limit";
    return out;
}
} // namespace ihb
#endif
