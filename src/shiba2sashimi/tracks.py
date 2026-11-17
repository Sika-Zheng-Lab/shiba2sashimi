"""Sample selection and arithmetic group summaries, independent of rendering."""
import numpy as np


def select_samples(experiment, samples=None, groups=None):
    def names(value):
        result = [item.strip() for item in value.split(",")]
        if any(not item for item in result) or len(set(result)) != len(result):
            raise ValueError("Sample/group selections must be nonempty and unique")
        return result

    if groups:
        order = names(groups)
        missing = set(order) - {info["group"] for info in experiment.values()}
        if missing:
            raise ValueError(f"Groups not found in experiment table: {', '.join(sorted(missing))}")
        return [sample for group in order for sample, info in experiment.items() if info["group"] == group]
    if samples:
        order = names(samples)
        missing = set(order) - set(experiment)
        if missing:
            raise ValueError(f"Samples not found in experiment table: {', '.join(sorted(missing))}")
        return order
    return list(experiment)


def build_tracks(coverage, junctions, experiment, order, psi=None, group_average=False):
    """Mean raw depths/counts include zero observations; PSI omits missing values.

    PSI is the unweighted mean of valid sample PSI, not Shiba's pooled or
    regression-predicted group PSI. Absent junction IDs count as zero; absent
    sample columns are an input error and must be checked by the caller.
    """
    members = {}
    for sample in order:
        key = experiment[sample]["group"] if group_average else sample
        members.setdefault(key, []).append(sample)
    result = []
    for label, samples in members.items():
        all_junctions = dict.fromkeys(j for sample in samples for j in junctions.get(sample, {}))
        values = []
        for sample in samples:
            try:
                value = float((psi or {}).get(sample, "NA"))
            except (ValueError, TypeError):
                continue
            if np.isfinite(value):
                values.append(value)
        result.append({
            "label": label, "group": experiment[samples[0]]["group"], "n": len(samples),
            "coverage": np.mean([coverage[sample] for sample in samples], axis=0),
            "junctions": {j: sum(junctions.get(s, {}).get(j, 0) for s in samples) / len(samples) for j in all_junctions},
            "psi": float(np.mean(values)) if values else None,
            "psi_n": len(values), "has_psi": psi is not None,
        })
    return result
