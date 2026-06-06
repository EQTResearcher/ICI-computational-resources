from ici import ICIParameters, compute_ici

params = ICIParameters.from_domain_defaults(
    domain="biological_lab",
    D=2.0e3,
    C=1.0e8,
    S=1.0e5,
    F=20,
    W=1.0e-3,
    M=1.0e2,
    system_name="E. coli K-12",
)

result = compute_ici(params)
print(result.summary())
