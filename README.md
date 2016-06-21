# mediation_OAI

####rand_int : implements the random interventional analogue to direct and indirect natural effects <br>
####ipw_med : implements calculations of controlled direct effect with IPTW, optionally setting M to different values or distributions. <br>
####med_reg : implements regression-based calculation of CDE, NDE, and NIE<br><br>

Everything assumes <br>
1. No X -> Y unmeasured confounder <br>
2. No M -> Y unmeasured confounder <br><br>

Addtionally<br>
3. med_reg assumes no exposure-induced M -> Y confounder <br>
4. NDE and their randomized analogue assume no X -> M confounder <br>


###Reference.
VanderWeele, T. (2015). Explanation in causal inference: methods for mediation and interaction. Oxford University Press.
