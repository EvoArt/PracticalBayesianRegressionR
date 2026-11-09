# PracticalBayesianRegressionR

Bayesian regression in R with the syntax of [brms](https://paul-buerkner.github.io/brms/),
fitted by [PracticalBayesianRegression.jl](https://github.com/EvoArt/PracticalBayesianRegression.jl).
The Julia package builds each model's log density and its gradient in closed form, straight
from the formula, and samples them with its own NUTS. You do not need Stan or a C++ toolchain.

```r
library(PracticalBayesianRegressionR)

fit <- brm(count ~ zAge + zBase * Trt + (1 | patient), data = epilepsy,
           family = poisson(), prior = prior(normal(0, 1), class = b))
summary(fit)
fixef(fit)
pp <- posterior_predict(fit)
loo(fit)
```

`bf()`, `lf()`, `set_rescor()`, `prior()`, `set_prior()`, `brmsfamily()` and the family
constructors (`student()`, `negbinomial()`, `zero_inflated_poisson()`, `cumulative()`, ...)
work as in brms. Fits convert to [posterior](https://mc-stan.org/posterior/) draws
(`as_draws_df(fit)`), so bayesplot and friends work too.

## Two back ends

**Julia** (`backend = "julia"`). A Julia session runs inside R through JuliaCall. It starts the
first time something needs it, and you can configure it with `pbr_setup()`. It covers every
model the Julia package supports.

**Compiled models** (`backend = "binary"`). A model is compiled once into a standalone
executable, which takes a few minutes and needs Julia and JuliaC. The executable is kept in a
model store. After that, `brm()` refits the model to any data with the right columns
**without starting Julia**: it writes the data as CSV, runs the executable and reads the
draws back. Given the same `seed`, the draws are identical to the Julia back end's.

```r
m <- compile_model(y ~ x + (1 | g), data = d, family = bernoulli())   # once
fit <- brm(m, data = new_data, seed = 1)                                # any session, no Julia
fit2 <- update(fit, newdata = more_data)                                # reruns the executable

compiled_models()                 # what is in the store

# by name, as a function
compile_model(y ~ x + (1 | g), data = d, family = bernoulli(), name = "mymod")
mymod <- load_model("mymod")      # or load_model("path/to/mymod/bin/mymod.exe")
fit <- mymod(new_data, chains = 4, seed = 1)
```

A model's directory is self-contained (executable, libraries, `model.toml`): copy the whole
directory to move a model to another machine with the same OS, then load it by its path.

The default, `backend = "auto"`, uses a compiled model when the store already has one for this
formula, family, priors and column types, and Julia otherwise. So once a model is compiled,
plain `brm(...)` calls with that model are fast in every later session.

Not every model can be compiled yet. Ordinal and categorical families, smooths, GPs,
censoring, autocorrelation and similar features are refused by name. Those models still fit
with the Julia back end.

## The model store is shared with Julia

Compiled models live in `~/.pbr/models`; set `PBR_MODEL_STORE` to use another directory. The
store is shared with the Julia package. A model compiled from R appears in Julia's
`compiled_models()` and runs with `brm(compiled_model("name"), data)`, and a model compiled
from Julia appears here in the same way. Each model is a directory with a `model.toml` that
describes it, plus the executable bundle.

## Set-up

- Julia 1.12 or later. JuliaCall finds it on the `PATH` or through `JULIA_HOME`.
- `pbr_setup(dev = "~/.julia/dev/PracticalBayesGLM")` uses a local checkout of the Julia
  package; without `dev`, it installs the GitHub version into a project of its own
  (`pbr_julia_project()`, or set `PBR_JULIA_PROJECT`).
- Compiling also needs a C compiler and JuliaC: either the `juliac` app on the `PATH`
  (`pkg> app add JuliaC`), or a command line in `options(pbr.juliac = ...)` or the
  `PBR_JULIAC` environment variable, for example
  `julia --project=/path/juliac-env -e "using JuliaC; JuliaC.main(ARGS)" --`.

## Differences from brms

- Functions inside formulas (`log(x)`, `I(x^2)`) are evaluated by Julia, so use functions
  Julia has under the same name. Variable names may contain only letters, digits and `_`.
- Text and factor columns are coded with their levels sorted, as `factor()` does by default.
  A factor whose levels are in a different order triggers a warning, because its reference
  level would differ from what brms uses.
- `brmsfamily("gamma")` without a link uses the log link. `Gamma()` from stats carries its
  inverse link explicitly, as in brms.
- `control` takes `adapt_delta`, `max_treedepth` and `metric` (`"diag"` or `"dense"`).

## Tests

```sh
Rscript -e 'testthat::test_local()'                    # no Julia needed
PBR_TEST_JULIA=1 Rscript -e 'testthat::test_local()'   # + Julia back end
PBR_TEST_JULIA=1 PBR_TEST_JULIAC=1 ...                 # + compiling (minutes)
```
