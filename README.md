# mediation_OAI

rand_int : implements the random interventional analogue to direct and indirect natural effects
ipw_med : implements calculations of controlled direct effect with IPTW, optionally setting M to different values or distributions.
med_reg : implements regression-based calculation of CDE, NDE, and NIE

Everything assumes 
(1) No X -> Y unmeasured confounder
(2) No M -> Y unmeasured confounder

Addtionally
(3) med_reg assumes no exposure-induced M -> Y confounder
(4) NDE and their randomized analogue assume no X -> M confounder
