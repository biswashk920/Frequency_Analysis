"""floodfreq: flood frequency analysis with Gumbel, GEV, Lognormal and Log-Pearson III."""
from .distributions import fit_gumbel, fit_gev, fit_lognormal, fit_lp3, lmoments
from .fitting import FITTERS, RETURN_PERIODS, fit_all, quantiles, plotting_positions, gof
from .bootstrap import bootstrap_ci

__version__ = "0.1.0"
