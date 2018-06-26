## Bayesian estimation of chronic disease epidemiology from incomplete data 


## Theoretical disease model 

We represent a disease as a continuous-time, multi-state process with three states: 

1. disease-free
2. disease
3. dead from the disease.

If we assume that mortality from other causes is independent of disease status, we can ignore deaths from other causes.   Assume also that recovery from the disease is not possible, but the framework is easily extensible to handle recurrence. The disease process is then fully defined by the 

* disease incidence $i(t)$, and the 

* case fatality rate $f(t)$,

which are both assumed to depend on age $t$.  Assume further that these rates are constant within integer years of age $a$, so they can be written $i_a$, $f_a$.

From these we can determine the _transition probability matrix_ $P_a$, whose $r,s$ entry is the probability that a person is in state $s$ at age $a+1$ given they are in state $r$ at age $a$. The matrix $P_a$ is defined in terms of $i_a$ and $f_a$ as a differential equation, whose solution is written out explicitly in the DisMod II paper. 

Further, let $p_a$ be the proportion of individuals in a birth cohort who are in each state at age $a$.  This is a row vector with three elements, one for each state $r$.  The transition probability matrix can then be used to obtain the corresponding proportions at age $a+1$, as 

$$ p_{a+1} = p_a P_a $$ 

The prevalence of disease (among people who are alive) at each age $a$ is then obtained as $pr_a = p_{a2} / (p_{a1} + p_{a2})$. 

The disease-specific mortality rate at age $a$, or the probability that a person alive at age $a$ dies from the disease before age $a+1$, can also be expressed in terms of the transition probabilities and prevalence, as 

$$ dm_a = P_a23 pr_a + P_a13 (1 - pr_a) $$


## Bayesian approach to estimating a model from data

Data are observed which give information about some, but not all, of the parameters in the theoretical disease model.   The form of the data available may be different in each application.  We then wish to estimate any remaining unknown parameters.   The Bayesian approach to estimating these unknowns can be described as four steps:

1. write down a theoretical model for underlying disease progression (as done above)

2. write down a statistical model for the observed data given the parameters of the underlying disease model.  An example of this is illustrated below.

3. write down prior distributions for the unknowns.  These may express prior ignorance, as in the example below.

4. compute the (unique) posterior distribution of the unknown parameters in the joint model, given the observed data.  Software discussed below. 

This approach is used in the DisMod-MR software, as explained in the book by Flaxman et al, however the software itself is undocumented and not fully published.  The older (published) DisMod II used an optimisation approach to estimate parameters, which is a rough less statistically principled approximation.  Other advantages of the Bayesian method include 

* straightforward to build in multiple sources of direct/indirect data.  This is enabled by the computational methods and available software for Bayesian modelling (illustrated below).   This allows the approach to generalise to settings with different forms of data available.  In contrast, DisMod II only allows limited forms of data. 

* automatic quantification of uncertainty, given the information supplied (as illustrated below)




## Example: IHD in Greater London 

The following data are given in the example spreadsheet 

* incidence of IHD by age

* IHD-specific mortality by age

and we wish to estimate case fatality by age, given these inputs. 

Note that these "data" are themselves estimates, based on underlying data from a finite population, that are not given.  In the Bayesian approach we can acknowledge that the incidence, mortality and/or case fatality are unknowns.   We illustrate this approach.  For simplicity, suppose that 

* the incidence $i_a$ is known and fixed at the quantities supplied in the spreadsheet, but 

* the disease-specific mortalities $dm_a$ are unknown, and have been estimated in terms of counts of individuals $n_a$ alive at age $a$ and the corresponding number $r_a$ who die of IHD before age $a+1$.

In this example we can roughly recreate $n_a$ and $r_a$ as follows.  We are given population counts for five-year age groups, thus we obtain approximate one-year population counts n_a as the dividing these values by 5.  The number of IHD deaths r_a is then approximately reconstructed by multiplying n_a by the IHD-specific mortality rates supplied in the spreadsheet. 

The statistical model is r_a ~ Binomial(n_a, m_a) 

where m_a is the probability of death from the disease 

m_a = P[1,3]*(1 - pr[a]) + P[2,3]*pr[a]

pr[a] computed given the P and the 
