"""Internal 21-mer supplement to exact junction evidence ("A" in names).

* ``acquisition``    internal-only safe 21-mer codebook projected on the
  msBWT; one-mate family-resolved pairs get the same finite host / reference
  state as the terminal branch (reference completeness judged on the pair's
  observed internal probes) and yield internal exact-junction witnesses;
* ``observability``  under *fixed* terminal proposals, per 3/2 fold
  rotation, exact events confirmed with existing plus internal A evidence
  but not with existing A evidence alone; same-family 500-bp locus
  components classified against prior validated full-memory events and
  cross-family local competition.

Internal 21-mers never create or widen proposals or confirm an event by
themselves.  They only make "A" observable.
"""
