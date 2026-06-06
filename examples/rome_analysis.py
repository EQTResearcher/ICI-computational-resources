from ici import ICIParameters, compute_ici_timeseries
from ici.visualization import plot_rt_timeseries
import numpy as np


def main():
    time_labels = [
        "Peak (100 CE)",
        "Third-century crisis (235 CE)",
        "Diocletian reform (285 CE)",
        "Late empire (400 CE)",
        "Western fall (476 CE)",
    ]
    params_series = [
        ICIParameters.from_domain_defaults("historical_empire", D=50, C=5e7, S=1e6, F=20, W=4.8e-13, M=1e3, system_name="Rome peak"),
        ICIParameters.from_domain_defaults("historical_empire", D=65, C=5e7, S=8e5, F=12, W=2.4e-13, M=6e2, system_name="Rome crisis"),
        ICIParameters.from_domain_defaults("historical_empire", D=55, C=5e7, S=9e5, F=18, W=4.0e-13, M=9e2, system_name="Diocletian reform"),
        ICIParameters.from_domain_defaults("historical_empire", D=80, C=3e7, S=4e5, F=6, W=9.5e-14, M=3e2, system_name="Late empire"),
        ICIParameters.from_domain_defaults("historical_empire", D=90, C=1.5e7, S=2e5, F=3, W=3.2e-14, M=1e2, system_name="Western fall"),
    ]
    results = compute_ici_timeseries(params_series, reference_index=0)
    for label, res in zip(time_labels, results):
        print(f"{label:30s} | ICI={res.ICI:6.1f} | R(t)={res.R_t:+.2f} | alert={res.alert_level.upper()}")

    fig = plot_rt_timeseries(
        time_points=list(range(len(time_labels))),
        Rt_values=np.array([r.R_t for r in results]),
        Rt_lower=np.array([r.R_t_lower for r in results]),
        Rt_upper=np.array([r.R_t_upper for r in results]),
        system_name="Roman Empire",
        events={1: "Crisis", 2: "Reform", 4: "Fall"},
    )
    fig.savefig("rome_rt_trajectory.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
