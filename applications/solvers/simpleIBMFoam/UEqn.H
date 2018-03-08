// Solve the Momentum equation
tmp<fvVectorMatrix> tUEqn
(
    new fvVectorMatrix
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    )
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);
solve(UEqn + fvc::grad(p));
phi = fvc::flux(U);

{
    UEqn -=
        particles.forcing
        (
            fvc::div(phi, U)
          + fvc::grad(p)
          + (turbulence->divDevReff(U) & U)
        );
}

UEqn.relax();
fvOptions.constrain(UEqn);
