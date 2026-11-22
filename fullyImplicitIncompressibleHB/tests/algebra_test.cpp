#include "../solver/HBAlgebra.H"
#include <cassert>
#include <iostream>

static void near(double a, double b, double tolerance)
{
    if (std::abs(a-b) > tolerance)
    {
        std::cerr << "Expected " << a << " near " << b << std::endl;
        std::abort();
    }
}

int main()
{
    near(ihb::rms(ihb::Vec{1e200, 1e200})/1e200, 1.0, 1e-15);
    near(ihb::rms(ihb::Vec{1e-200, 1e-200})/1e-200, 1.0, 1e-15);
    const double pi = std::acos(-1.0);
    for (int h = 0; h <= 8; ++h)
    {
        const double omega = 1.7;
        const std::vector<ihb::Vec> d = ihb::differentiation(h, omega);
        const int n = 2*h+1;
        for (int i = 0; i < n; ++i)
        {
            near(ihb::dot(d[i], ihb::Vec(n, 1.0)), 0.0, 1e-12);
            for (int j = 0; j < n; ++j) near(d[i][j], -d[j][i], 1e-15);
            for (int k = 1; k <= h; ++k)
            {
                ihb::Vec sine(n), cosine(n);
                for (int j = 0; j < n; ++j)
                {
                    sine[j] = std::sin(2*pi*k*j/n);
                    cosine[j] = std::cos(2*pi*k*j/n);
                }
                near(ihb::dot(d[i], sine), omega*k*cosine[i], 1e-11);
                near(ihb::dot(d[i], cosine), -omega*k*sine[i], 1e-11);
            }
        }
    }
    bool threw = false;
    try { ihb::differentiation(-1, 1.0); }
    catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    // Indefinite saddle point system, with exact solution (1, -2, 3).
    ihb::Vec x(3, 0.0);
    const ihb::Residual saddle = [](const ihb::Vec& a)
    {
        return ihb::Vec{2*a[0]+a[2]-5, 3*a[1]+a[2]+3, a[0]+a[1]+1};
    };
    ihb::Settings s;
    s.absTol = 1e-10;
    s.relTol = 1e-10;
    s.linearRelTol = 1e-6;
    ihb::Result result = ihb::newton(x, saddle, s);
    assert(result.converged);
    near(x[0], 1.0, 1e-8); near(x[1], -2.0, 1e-8); near(x[2], 3.0, 1e-8);

    // Large first Newton step requires backtracking for x^2=2.
    x = ihb::Vec{0.1, 0.5};
    double smallestAlpha = 1.0;
    result = ihb::newton(x, [](const ihb::Vec& a)
    {
        return ihb::Vec{a[0]*a[0]-2.0, std::exp(a[1])-3.0};
    }, s, [&](int i, double, double alpha)
    {
        if (i) smallestAlpha = std::min(smallestAlpha, alpha);
    });
    assert(result.converged && smallestAlpha < 1.0);
    near(x[0], std::sqrt(2.0), 1e-8); near(x[1], std::log(3.0), 1e-8);

    // Restarting actually works (more iterations than Krylov basis size).
    ihb::Vec b(12, 1.0);
    const ihb::LinearResult restarted = ihb::gmres([](const ihb::Vec& a)
    {
        ihb::Vec y(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) y[i] = (1.0+0.1*i)*a[i];
        return y;
    }, b, 3, 200, 1e-11, 1e-12);
    assert(restarted.converged && restarted.iterations > 3);
    for (std::size_t i = 0; i < b.size(); ++i)
        near(restarted.x[i], 1.0/(1.0+0.1*i), 1e-10);

    x = ihb::Vec{0.1, 0.5};
    ihb::Settings limited = s;
    limited.maxNewton = 0;
    result = ihb::newton(x, [](const ihb::Vec& a)
    { return ihb::Vec{a[0]*a[0]-2.0, a[1]-1.0}; }, limited);
    assert(!result.converged && result.reason == "Newton iteration limit");

    // Singular inconsistent problem cannot be marked converged.
    x = ihb::Vec{1.0};
    result = ihb::newton(x, [](const ihb::Vec&)
    { return ihb::Vec{1.0}; }, s);
    assert(!result.converged);

    result = ihb::newton(x, [](const ihb::Vec&)
    { return ihb::Vec{std::numeric_limits<double>::quiet_NaN()}; }, s);
    assert(!result.converged && result.evaluations == 1);

    // Reject dimension mismatch instead of indexing beyond a buffer.
    result = ihb::newton(x, [](const ihb::Vec&)
    { return ihb::Vec{1.0, 2.0}; }, s);
    assert(!result.converged);
    std::cout << "PASS: Fourier constant/skew/sine/cosine; saddle point; "
              << "nonlinear/backtracking; restarted GMRES; failure paths\n";
}
