"""
Architecture extraction notes for AI ICI cases.

Recommended extraction components:
- computational graph parsing for D;
- effective model/data capacity estimation for C;
- measured hardware throughput or roofline-corrected throughput for S;
- directed cycle / recurrent / attention-feedback analysis for F;
- state-refresh or online-update frequency for W;
- effective memory complexity after redundancy and Hessian-condition correction for M.
"""

RECOMMENDED_ARTIFACTS = [
    "architecture graph",
    "parameter-count and active-parameter estimate",
    "effective data/model capacity estimate",
    "measured throughput or roofline-corrected throughput",
    "feedback-loop / recurrence / attention-integration estimate",
    "state-refresh / online-adaptation estimate",
    "memory-complexity estimate",
]
