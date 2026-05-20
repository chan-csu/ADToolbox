# ADM Models

ADToolbox exposes two anaerobic digestion models:

| CLI command | Model | Parameter folder | Stoichiometric builder | ODE system |
| --- | --- | --- | --- | --- |
| `adtoolbox ADM adm1` | ADM1 | `adm1` | `adm.build_adm1_stoichiometric_matrix` | `adm.adm1_ode_sys` |
| `adtoolbox ADM e-adm` | e-ADM | `e_adm` | `adm.build_e_adm_stoichiometric_matrix` | `adm.e_adm_ode_sys` |

The CLI and public API use only the names ADM1 and e-ADM. The preferred parameter format is one `models.json` file containing all models, keyed by model name.

## Model Inputs

The canonical model input is a single JSON file with one top-level key per model:

```json
{
  "adm1": {
    "species": [],
    "reactions": [],
    "model_parameters": {},
    "base_parameters": {},
    "initial_conditions": {},
    "inlet_conditions": {}
  },
  "e_adm": {
    "species": [],
    "reactions": [],
    "model_parameters": {},
    "base_parameters": {},
    "initial_conditions": {},
    "inlet_conditions": {}
  }
}
```

The first layer is always the model key. The second layer is the complete model payload used to instantiate `adm.Model`.

| Top-level key | Public model | CLI command |
| --- | --- | --- |
| `adm1` | ADM1 | `adtoolbox ADM adm1 --models-json models.json` |
| `e_adm` | e-ADM | `adtoolbox ADM e-adm --models-json models.json` |

The older six-file layout is still supported as a compatibility path, either by passing `--parameters-dir` or by passing each file explicitly.

| Input file | ADM1 name | e-ADM name | Description |
| --- | --- | --- | --- |
| Model parameters | `adm1_model_parameters.json` | `e_adm_model_parameters.json` | Kinetic, equilibrium, inhibition, gas transfer, and yield parameters used in reaction rates and acid/base calculations. |
| Base parameters | `adm1_base_parameters.json` | `e_adm_base_parameters.json` | Reactor-level constants such as liquid volume, gas volume, influent flow, gas transfer constants, pressure, temperature, and gas constants. |
| Initial conditions | `adm1_initial_conditions.json` | `e_adm_initial_conditions.json` | Initial concentration of every model state. Keys must match the species list. |
| Inlet conditions | `adm1_inlet_conditions.json` | `e_adm_inlet_conditions.json` | Influent concentration of every model state, using the same species names with `_in` appended. |
| Reactions | `adm1_reactions.json` | `e_adm_reactions.json` | Ordered reaction names. This order defines the columns of the stoichiometric matrix and the reaction-rate vector. |
| Species | `adm1_species.json` | `e_adm_species.json` | Ordered state names. This order defines the rows of the stoichiometric matrix and the concentration vector. |

The CLI also accepts:

| Option | Applies to | Description |
| --- | --- | --- |
| `--models-json` | ADM1, e-ADM | Preferred input. JSON file containing all models keyed by model name. |
| `--parameters-dir` | ADM1, e-ADM | Directory containing the six JSON files for the selected model. |
| `--model-parameters`, `--base-parameters`, `--initial-conditions`, `--inlet-conditions`, `--reactions`, `--species` | ADM1, e-ADM | Override individual JSON file paths. |
| `--report` | ADM1, e-ADM | Output mode. Use `csv` to write a report or `dash` to open the interactive app. If omitted, the model is solved without creating an output view. |
| `--control-states` | e-ADM | JSON object of states to hold constant. By default, e-ADM holds `S_H_ion` at `10^-6.5`. |

### ADM1 Input Contract

| Input object | Required shape | Required naming rule | How it is used |
| --- | --- | --- | --- |
| `species` | list of 38 strings | Names must match every key in `initial_conditions`; each inlet key must be `species_name + "_in"`. | Defines the order of the state vector `c`. |
| `reactions` | list of 28 strings | Names must match the reaction names used in `build_adm1_stoichiometric_matrix` and `adm1_ode_sys`. | Defines the order of the reaction vector `v` and the columns of `S`. |
| `initial_conditions` | object with 38 numeric values | Keys are exactly the ADM1 species names. | Converted to `y0`, the initial state passed to `solve_ivp`. |
| `inlet_conditions` | object with 38 numeric values | Keys are exactly the ADM1 species names with `_in` appended. | Converted to `c_in`, used in the liquid and gas flow terms. |
| `base_parameters` | object with reactor constants | Must include `q_in`, `V_liq`, `V_gas`, `R`, `T_op`, `T_base`, `P_atm`, and `K_W` or equivalent water-equilibrium value used by the model. | Controls hydraulic dilution, gas headspace scaling, and gas partial pressures. |
| `model_parameters` | object with kinetic, yield, fraction, acid/base, inhibition, and gas-transfer constants | Must include all keys referenced by ADM1 rates and stoichiometry, such as `k_dis`, `k_hyd_*`, `k_m_*`, `K_S_*`, `Y_*`, `f_*`, `C_*`, `N_*`, `K_a_*`, `K_pH_*`, `K_I_*`, `k_L_a`, `K_H_*`, and `k_p`. | Controls reaction rates, stoichiometric coefficients, inhibition factors, acid/base rates, and gas transfer. |

### e-ADM Input Contract

| Input object | Required shape | Required naming rule | How it is used |
| --- | --- | --- | --- |
| `species` | list of 49 strings | Names must match every key in `initial_conditions`; each inlet key must be `species_name + "_in"`. | Defines the order of the state vector `c`. |
| `reactions` | list of 45 strings | Names must match the reaction names used in `build_e_adm_stoichiometric_matrix` and `e_adm_ode_sys`. | Defines the order of the reaction vector `v` and the columns of `S`. |
| `initial_conditions` | object with 49 numeric values | Keys are exactly the e-ADM species names. | Converted to `y0`, after applying any `control_state` overrides. |
| `inlet_conditions` | object with 49 numeric values | Keys are exactly the e-ADM species names with `_in` appended. | Converted to `c_in`, used in the liquid and gas flow terms. |
| `base_parameters` | object with reactor constants | Must include `q_in`, `V_liq`, `V_gas`, `R`, `T_op`, `T_base`, `P_atm`, and `K_W`. | Controls hydraulic dilution, gas headspace scaling, and gas partial pressures. |
| `model_parameters` | object with kinetic, yield, fraction, acid/base, inhibition, and gas-transfer constants | Must include all keys referenced by e-ADM rates and stoichiometry, including `Y_su`, `Y_aa`, `Y_fa`, `Y_ac_et`, `Y_ac_lac`, `Y_pro_et`, `Y_pro_lac`, `Y_bu_et`, `Y_bu_lac`, `Y_va`, `Y_cap`, `Y_bu`, `Y_Me_ac`, `Y_Me_CO2`, `Y_ac_et_ox`, and `Y_pro_lac_ox`. | Controls reaction rates, stoichiometric coefficients, inhibition factors, acid/base rates, and gas transfer. |
| `control_state` | optional object | Keys must be valid species names. | Overrides initial conditions and forces selected derivatives to zero. The CLI sets `S_H_ion = 10^-6.5` by default for e-ADM. |

### Other Keyed Database Files

The same first-layer-key pattern should be used for other structured databases. For example, feeds can live in one `feeds.json` file:

```json
{
  "food_waste": {
    "name": "food_waste",
    "carbohydrates": 10,
    "proteins": 20,
    "lipids": 20,
    "si": 30,
    "xi": 50,
    "tss": 80
  },
  "manure": {
    "name": "manure",
    "carbohydrates": 5,
    "proteins": 15,
    "lipids": 5,
    "si": 20,
    "xi": 55,
    "tss": 70
  }
}
```

That keeps the rule consistent: one file per database type, one top-level key per named item.

## Model Outputs

| Output | Source | Description |
| --- | --- | --- |
| Concentration trajectories | `Model.solve_model(...)` | Returns a SciPy `OdeResult`. `sol.t` contains time points and `sol.y` contains species concentrations ordered exactly like `model.species`. |
| CSV report | `--report csv` or `Model.csv_report(...)` | Writes `<model.name>_Report.csv`. Rows are species and columns are simulated time points. |
| Interactive dashboard | `--report dash` or `Model.dash_app(...)` | Opens a Dash app with concentration plots and model tables. |
| Reaction fluxes | `model.info["Fluxes"]` | The most recent reaction-rate vector `v`, ordered exactly like `model.reactions`. |
| Stoichiometric matrix | `model.s` | Matrix `S` with shape `number_of_species x number_of_reactions`. Entry `S[i, j]` is the coefficient for species `i` in reaction `j`. |

ADM1 has 38 species and 28 reactions. e-ADM has 49 species and 45 reactions.

### Output Shapes

| Output | ADM1 shape | e-ADM shape | Notes |
| --- | --- | --- | --- |
| `sol.t` | `(n_time_points,)` | `(n_time_points,)` | Time is in days in the default examples. The CLI uses `np.linspace(0, 30, 10000)`. |
| `sol.y` | `(38, n_time_points)` | `(49, n_time_points)` | Rows follow `model.species`; columns follow `sol.t`. |
| `model.s` | `(38, 28)` | `(49, 45)` | Rows follow species; columns follow reactions. |
| `model.info["Fluxes"]` | `(28, 1)` | `(45, 1)` | Stores the most recent flux calculation from the ODE function. |
| CSV report | 38 rows plus header | 49 rows plus header | The first column is the species index. Remaining columns are simulated time points. |

The solver returns a fallback solution filled with `1e10` if integration fails. That is a failure sentinel, not a physical concentration.

## Mathematical Setup

The JSON files describe the state names, reaction names, and parameter values. The Python model converts them into a dynamic system by building a stoichiometric matrix and evaluating reaction rates at every time step.

| Step | Formula | Meaning |
| --- | --- | --- |
| Species vector | `c(t) = [c_1(t), ..., c_n(t)]^T` | Concentrations are ordered by the species JSON file. |
| Reaction-rate vector | `v(c, theta) = [v_1, ..., v_m]^T` | Fluxes are ordered by the reactions JSON file and are calculated from concentrations and parameters. |
| Stoichiometric matrix | `S in R^(n x m)` | Rows are species, columns are reactions. `S[i, j]` converts reaction `j` into the concentration change of species `i`. |
| Reaction contribution | `r(c, theta) = S v(c, theta)` | Net production or consumption of every species from biochemical, acid/base, and gas-transfer reactions. |
| Liquid flow contribution | `(q_in / V_liq) * (c_in - c)` | Continuous stirred tank dilution or influent loading for liquid-phase states. |
| Gas flow contribution | `(q_gas / V_gas) * (c_gas,in - c_gas)` | Headspace gas outflow for gas states. |
| Gas flow rate | `q_gas = max(0, k_p * (P_gas - P_atm))` | Gas leaves only when calculated gas pressure exceeds atmospheric pressure. |
| Final derivative | `dc/dt = S v(c, theta) + flow terms` | This is the vector returned to `scipy.integrate.solve_ivp`. |

The implemented derivative can be written as:

```text
c = current state vector
c_in = inlet state vector
S = stoichiometric matrix
v = reaction-rate vector

r = S @ v

for liquid states:
    dc_liq/dt = r_liq + (q_in / V_liq) * (c_in,liq - c_liq)

for gas states:
    dc_gas/dt = r_gas + (q_gas / V_gas) * (c_in,gas - c_gas)
```

The gas states are the final three species in both models:

```text
S_gas_h2, S_gas_ch4, S_gas_co2
```

The gas pressure and gas outflow equations are:

```text
p_gas_h2  = S_gas_h2  * R * T_op / 16
p_gas_ch4 = S_gas_ch4 * R * T_op / 64
p_gas_co2 = S_gas_co2 * R * T_op
p_gas_h2o = 0.0313 * exp(5290 * (1 / T_base - 1 / T_op))

P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o
q_gas = max(0, k_p * (P_gas - P_atm))
```

For any reaction `j`, the instantaneous contribution to species `i` is:

```text
dc_i/dt from reaction j = S[i, j] * v_j
```

This is why reaction and species names are used to build `S`: the equations are much easier to audit than when the matrix is filled with positional numbers.

## Reaction Rate Forms

| Reaction type | Formula pattern | Used for |
| --- | --- | --- |
| First-order conversion | `v_j = k * X` | Disintegration, hydrolysis, and decay reactions. |
| Single-substrate uptake | `v_j = k_m * S_a / (K_S + S_a) * X * I` | Sugar, amino acid, LCFA, acetate, hydrogen, ethanol, lactate, and VFA uptake reactions. |
| Competitive ADM1 C4 uptake | `v_j = k_m * S_a / (K_S + S_a) * X_c4 * S_a / (S_va + S_bu + eps) * I` | ADM1 valerate and butyrate uptake. |
| Two-substrate uptake | `v_j = k_m * S_a * S_b / (K_a * S_a + K_b * S_b + S_a * S_b + eps) * X * I` | e-ADM syntrophic reactions that require paired substrates. |
| Inhibition or limitation | `I = product(I_pH, I_H2, I_NH3, I_IN, ...)` | Multiplies uptake rates to suppress a reaction under limiting or inhibitory conditions. |
| Acid/base kinetic form | `v_AB = k_AB * (S_ion * (K_a + S_H) - K_a * S_total)` | Kinetic acid/base reactions when the state is not solved algebraically. |
| DAE acid/base update | `S_ion = K_a / (K_a + S_H) * S_total` | Algebraic speciation used when `model.switch == "DAE"`. |
| Gas transfer | `v_gas = k_La * (S_liq - H * p_gas)` | Transfer between dissolved gas states and headspace gas states. |

### Acid/Base Equations

ADM1 and e-ADM both calculate kinetic acid/base rates, then set selected acid/base derivatives to zero in DAE mode. The kinetic form is:

```text
v_AB,acid = k_A_B,acid * (S_ion * (K_a,acid + S_H) - K_a,acid * S_total)
```

For carbonate and inorganic nitrogen:

```text
v_AB,co2 = k_A_B_co2 * (S_hco3_ion * (K_a_co2 + S_H) - K_a_co2 * S_IC)
v_AB,IN  = k_A_B_IN  * (S_nh3      * (K_a_IN  + S_H) - K_a_IN  * S_IN)
```

In DAE mode, ionized species are updated algebraically. e-ADM explicitly updates:

```text
S_va_ion  = K_a_va  / (K_a_va  + S_H) * S_va
S_bu_ion  = K_a_bu  / (K_a_bu  + S_H) * S_bu
S_pro_ion = K_a_pro / (K_a_pro + S_H) * S_pro
S_cap_ion = K_a_cap / (K_a_cap + S_H) * S_cap
S_ac_ion  = K_a_ac  / (K_a_ac  + S_H) * S_ac
S_lac_ion = K_a_lac / (K_a_lac + S_H) * S_lac
S_hco3_ion = S_IC - S_co2
```

Hydrogen ion concentration is calculated from electroneutrality unless it is controlled:

```text
S_H = -phi / 2 + 0.5 * sqrt(phi^2 + 4 * K_w)
```

When `S_H_ion` is included in `control_state`, the controlled value is used and `dS_H_ion/dt = 0`.

### Inhibition Equations

The pH inhibition terms use Hill-type functions:

```text
I_pH_x = K_pH_x^n_x / (S_H^n_x + K_pH_x^n_x)
```

Nitrogen limitation is a Monod term:

```text
I_IN_lim = S_IN / (K_S_IN + S_IN + eps)
```

Hydrogen and ammonia inhibition terms use:

```text
I_h2_x = 1 / (1 + S_h2 / (K_I_h2_x + eps))
I_nh3  = 1 / (1 + S_nh3 / (K_I_nh3 + eps))
```

ADM1 uses the same conceptual forms, with its original `I_IN_lim = 1 / (1 + K_S_IN / S_IN)` expression.

## ADM1 Reactions

The ADM1 reactions follow the Batstone ADM1 structure: disintegration, hydrolysis, soluble-substrate uptake, biomass decay, acid/base speciation, and gas transfer.

### ADM1 Stoichiometric Coefficients

ADM1 uses yield coefficients `Y_*`, product fractions `f_*`, carbon contents `C_*`, and nitrogen contents `N_*` to fill `S`. Nitrogen limitation sets selected yields to zero before `S` is built.

For the main carbon-balance terms:

```text
Y_su = 0 if nitrogen_limited else Y_su
Y_aa = 0 if nitrogen_limited else Y_aa
Y_fa = 0 if nitrogen_limited else Y_fa

s_1  = -C_xc + f_sI_xc*C_sI + f_ch_xc*C_ch + f_pr_xc*C_pr + f_li_xc*C_li + f_xI_xc*C_xI
s_2  = -C_ch + C_su
s_3  = -C_pr + C_aa
s_4  = -C_li + (1 - f_fa_li)*C_su + f_fa_li*C_fa
s_5  = -C_su + (1 - Y_su)*(f_bu_su*C_bu + f_pro_su*C_pro + f_ac_su*C_ac) + Y_su*C_bac
s_6  = -C_aa + (1 - Y_aa)*(f_va_aa*C_va + f_bu_aa*C_bu + f_pro_aa*C_pro + f_ac_aa*C_ac) + Y_aa*C_bac
s_7  = -C_fa + (1 - Y_fa)*0.7*C_ac + Y_fa*C_bac
s_8  = -C_va + (1 - Y_c4)*0.54*C_pro + (1 - Y_c4)*0.31*C_ac + Y_c4*C_bac
s_9  = -C_bu + (1 - Y_c4)*0.8*C_ac + Y_c4*C_bac
s_10 = -C_pro + (1 - Y_pro)*0.57*C_ac + Y_pro*C_bac
s_11 = -C_ac + (1 - Y_ac)*C_ch4 + Y_ac*C_bac
s_12 = (1 - Y_h2)*C_ch4 + Y_h2*C_bac
s_13 = -C_bac + C_xc
```

These terms are inserted into the `S_IC` row with negative signs for conversion reactions and positive return through decay. For example:

```text
S[S_IC, Uptake of sugars] = -s_5
S[S_IC, Uptake of acetate] = -s_11
S[S_IC, Decay of Xsu] = -s_13
```

Nitrogen stoichiometry follows the same pattern:

```text
S[S_IN, Uptake of sugars] = -Y_su * N_bac
S[S_IN, Uptake of amino acids] = N_aa - Y_aa * N_bac
S[S_IN, Uptake of acetate] = -Y_ac * N_bac
S[S_IN, Decay of Xsu] = N_bac - N_xc
```

| Reaction | Rate expression | Main conversion represented in `S` |
| --- | --- | --- |
| `Disintegration` | `k_dis * X_xc` | Composite particulate matter `X_xc` is split into carbohydrates `X_ch`, proteins `X_pr`, lipids `X_li`, soluble inert `S_I`, particulate inert `X_I`, inorganic carbon `S_IC`, and inorganic nitrogen `S_IN`. |
| `Hydrolysis carbohydrates` | `k_hyd_ch * X_ch` | Particulate carbohydrates become soluble sugars `S_su`, with carbon correction through `S_IC`. |
| `Hydrolysis of proteins` | `k_hyd_pr * X_pr` | Particulate proteins become amino acids `S_aa`, with carbon correction through `S_IC`. |
| `Hydrolysis of lipids` | `k_hyd_li * X_li` | Lipids split into LCFA `S_fa` and sugars `S_su`, with carbon correction through `S_IC`. |
| `Uptake of sugars` | `k_m_su * S_su / (K_S_su + S_su) * X_su * I5` | Sugars are consumed to grow `X_su` and produce butyrate `S_bu`, propionate `S_pro`, acetate `S_ac`, hydrogen `S_h2`, inorganic carbon, and nitrogen demand. |
| `Uptake of amino acids` | `k_m_aa * S_aa / (K_S_aa + S_aa) * X_aa * I6` | Amino acids are consumed to grow `X_aa` and produce valerate `S_va`, butyrate, propionate, acetate, hydrogen, inorganic carbon, and inorganic nitrogen. |
| `Uptake of LCFA` | `k_m_fa * S_fa / (K_S_fa + S_fa) * X_fa * I7` | LCFA are consumed to grow `X_fa` and produce acetate, hydrogen, inorganic carbon, and nitrogen demand. |
| `Uptake of valerate` | `k_m_c4 * S_va / (K_S_c4 + S_va) * X_c4 * S_va / (S_va + S_bu + eps) * I8` | Valerate is consumed by `X_c4`, producing propionate, acetate, hydrogen, biomass, carbon, and nitrogen demand. |
| `Uptake of butyrate` | `k_m_c4 * S_bu / (K_S_c4 + S_bu) * X_c4 * S_bu / (S_bu + S_va + eps) * I9` | Butyrate is consumed by `X_c4`, producing acetate, hydrogen, biomass, carbon, and nitrogen demand. |
| `Uptake of propionate` | `k_m_pr * S_pro / (K_S_pro + S_pro) * X_pro * I10` | Propionate is consumed to grow `X_pro` and produce acetate, hydrogen, carbon, and nitrogen demand. |
| `Uptake of acetate` | `k_m_ac * S_ac / (K_S_ac + S_ac) * X_ac * I11` | Acetate is consumed by acetoclastic methanogens, producing methane `S_ch4`, `X_ac`, carbon correction, and nitrogen demand. |
| `Uptake of Hydrogen` | `k_m_h2 * S_h2 / (K_S_h2 + S_h2) * X_h2 * I12` | Hydrogen is consumed by hydrogenotrophic methanogens, producing methane `S_ch4`, `X_h2`, carbon correction, and nitrogen demand. |
| `Decay of Xsu` | `k_dec_X_su * X_su` | Sugar degraders decay back to composite particulate matter `X_xc`, with carbon and nitrogen returned through `S_IC` and `S_IN`. |
| `Decay of Xaa` | `k_dec_X_aa * X_aa` | Amino acid degraders decay to `X_xc`, with carbon and nitrogen returned. |
| `Decay of Xfa` | `k_dec_X_fa * X_fa` | LCFA degraders decay to `X_xc`, with carbon and nitrogen returned. |
| `Decay of Xc4` | `k_dec_X_c4 * X_c4` | C4 degraders decay to `X_xc`, with carbon and nitrogen returned. |
| `Decay of Xpro` | `k_dec_X_pro * X_pro` | Propionate degraders decay to `X_xc`, with carbon and nitrogen returned. |
| `Decay of Xac` | `k_dec_X_ac * X_ac` | Acetoclastic methanogens decay to `X_xc`, with carbon and nitrogen returned. |
| `Decay of Xh2` | `k_dec_X_h2 * X_h2` | Hydrogenotrophic methanogens decay to `X_xc`, with carbon and nitrogen returned. |
| `Acid Base Equilibrium (Va)` | `k_A_B_va * (S_va_ion * (K_a_va + S_H) - K_a_va * S_va)` | Valerate speciation between total valerate `S_va` and ionized valerate `S_va_ion`. |
| `Acid Base Equilibrium (Bu)` | `k_A_B_bu * (S_bu_ion * (K_a_bu + S_H) - K_a_bu * S_bu)` | Butyrate speciation between `S_bu` and `S_bu_ion`. |
| `Acid Base Equilibrium (Pro)` | `k_A_B_pro * (S_pro_ion * (K_a_pro + S_H) - K_a_pro * S_pro)` | Propionate speciation between `S_pro` and `S_pro_ion`. |
| `Acid Base Equilibrium (Ac)` | `k_A_B_ac * (S_ac_ion * (K_a_ac + S_H) - K_a_ac * S_ac)` | Acetate speciation between `S_ac` and `S_ac_ion`. |
| `Acid Base Equilibrium (CO2)` | `k_A_B_co2 * (S_hco3_ion * (K_a_co2 + S_H) - K_a_co2 * S_IC)` | Carbonate speciation through bicarbonate `S_hco3_ion` and inorganic carbon `S_IC`. |
| `Acid Base Equilibrium (In)` | `k_A_B_IN * (S_nh3 * (K_a_IN + S_H) - K_a_IN * S_IN)` | Nitrogen speciation between ammonia `S_nh3`, ammonium `S_nh4_ion`, and total inorganic nitrogen `S_IN`. |
| `Gas Transfer H2` | `k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2)` | Dissolved hydrogen leaves the liquid and enters `S_gas_h2` with `V_liq / V_gas` scaling. |
| `Gas Transfer CH4` | `k_L_a * (S_ch4 - 64 * K_H_ch4 * p_gas_ch4)` | Dissolved methane leaves the liquid and enters `S_gas_ch4` with `V_liq / V_gas` scaling. |
| `Gas Transfer CO2` | `k_L_a * (S_co2 - K_H_co2 * p_gas_co2)` | Dissolved carbon dioxide leaves the liquid and enters `S_gas_co2` with `V_liq / V_gas` scaling. |

### ADM1 Inhibition Factors

| Factor | Formula pattern | Applied to |
| --- | --- | --- |
| `I5` | `I_pH_aa * I_IN_lim` | Sugar uptake. |
| `I6` | `I5` | Amino acid uptake. |
| `I7` | `I_pH_aa * I_IN_lim * I_h2_fa` | LCFA uptake. |
| `I8`, `I9` | `I_pH_aa * I_IN_lim * I_h2_c4` | Valerate and butyrate uptake. |
| `I10` | `I_pH_aa * I_IN_lim * I_h2_pro` | Propionate uptake. |
| `I11` | `I_pH_ac * I_IN_lim * I_nh3` | Acetate uptake. |
| `I12` | `I_pH_h2 * I_IN_lim` | Hydrogen uptake. |

## e-ADM Reactions

e-ADM keeps the same dynamic structure but expands the soluble fermentation network. It separates ethanol-linked and lactate-linked routes, adds caproate, lactate, ethanol, explicit TSS/TDS disintegration, and two methanogenesis routes.

### e-ADM Stoichiometric Coefficients

e-ADM uses the same matrix logic as ADM1, but most soluble conversions are written as explicit product fractions plus a carbon-balance residual. Selected yields are set to zero under nitrogen limitation.

For sugar, amino acid, and LCFA uptake:

```text
Y_su = 0 if nitrogen_limited else Y_su
Y_aa = 0 if nitrogen_limited else Y_aa
Y_fa = 0 if nitrogen_limited else Y_fa

f_ac_su = 1 - f_pro_su - f_et_su - f_lac_su
f_ac_aa = 1 - f_pro_aa - f_et_aa - f_lac_aa
f_ac_fa = 1 - f_pro_fa - f_et_fa - f_lac_fa

f_IC_su = -(-C_su + (1-Y_su)*f_pro_su*C_pro + (1-Y_su)*f_et_su*C_et
              + (1-Y_su)*f_lac_su*C_lac + (1-Y_su)*f_ac_su*C_ac + Y_su*C_bac)

f_IC_aa = -(-C_aa + (1-Y_aa)*f_pro_aa*C_pro + (1-Y_aa)*f_et_aa*C_et
              + (1-Y_aa)*f_lac_aa*C_lac + (1-Y_aa)*f_ac_aa*C_ac + Y_aa*C_bac)

f_IC_fa = -(-C_fa + (1-Y_fa)*f_pro_fa*C_pro + (1-Y_fa)*f_et_fa*C_et
              + (1-Y_fa)*f_lac_fa*C_lac + (1-Y_fa)*f_ac_fa*C_ac + Y_fa*C_bac)
```

For chain elongation and VFA conversion, e-ADM uses paired substrate routes:

```text
f_IC_ac_et  = -(-C_ac  + f_et_ac*C_et  + (1-f_et_ac-Y_ac_et)*f_bu_ac*C_bu  + Y_ac_et*C_bac)
f_IC_ac_lac = -(-C_ac  + f_lac_ac*C_lac + (1-f_lac_ac-Y_ac_lac)*f_bu_ac*C_bu + Y_ac_lac*C_bac)

f_IC_pro_et  = -(-C_pro + f_et_pro*C_et  + (1-f_et_pro-Y_pro_et)*f_va_pro*C_va  + Y_pro_et*C_bac)
f_IC_pro_lac = -(-C_pro + f_lac_pro*C_lac + (1-f_lac_pro-Y_pro_lac)*f_va_pro*C_va + Y_pro_lac*C_bac)

f_IC_bu_et  = -(-C_bu + f_et_bu*C_et  + (1-f_et_bu-Y_bu_et)*f_cap_bu*C_cap  + Y_bu_et*C_bac)
f_IC_bu_lac = -(-C_bu + f_lac_bu*C_lac + (1-f_lac_bu-Y_bu_lac)*f_cap_bu*C_cap + Y_bu_lac*C_bac)
```

Methanogenesis and oxidation use:

```text
f_IC_Me_ach2 = -(f_ac_h2*C_ac + (1 - Y_Me_ac)*C_ch4 + Y_Me_ac*C_bac)
f_IC_et_ox   = -(-C_et  + (1-Y_ac_et_ox)*C_bac + Y_ac_et_ox*C_ac)
f_IC_lac_ox  = -(-C_lac + (1-Y_pro_lac_ox)*C_bac + Y_pro_lac_ox*C_pro)
```

The stoichiometric matrix then combines these residuals with explicit product coefficients. For example:

```text
S[S_su, Uptake of sugars] = -1
S[S_pro, Uptake of sugars] = (1 - Y_su) * f_pro_su
S[S_et, Uptake of sugars] = (1 - Y_su) * f_et_su
S[S_lac, Uptake of sugars] = (1 - Y_su) * f_lac_su
S[S_ac, Uptake of sugars] = (1 - Y_su) * f_ac_su
S[S_IN, Uptake of sugars] = -Y_su * N_bac
S[S_IC, Uptake of sugars] = f_IC_su
S[X_su, Uptake of sugars] = Y_su
```

| Reaction | Rate expression | Main conversion represented in `S` |
| --- | --- | --- |
| `TSS_Disintegration` | `k_dis_TSS * TSS` | Total suspended solids are split into `X_ch`, `X_pr`, `X_li`, and inert particulate biomass `X_I` using feed fractions. |
| `TDS_Disintegration` | `k_dis_TDS * TDS` | Total dissolved solids are split into `X_ch`, `X_pr`, `X_li`, and soluble inert `S_I` using feed fractions. |
| `Hydrolysis carbohydrates` | `k_hyd_ch * X_ch` | Particulate carbohydrate is hydrolyzed to soluble sugar `S_su`. |
| `Hydrolysis proteins` | `k_hyd_pr * X_pr` | Particulate protein is hydrolyzed to amino acids `S_aa`. |
| `Hydrolysis lipids` | `k_hyd_li * X_li` | Particulate lipids are hydrolyzed to LCFA `S_fa`. |
| `Uptake of sugars` | `k_m_su * S_su / (K_S_su + S_su) * X_su * I5` | Sugars are consumed to grow `X_su` and form propionate, ethanol, lactate, acetate, carbon, and nitrogen demand. |
| `Uptake of amino acids` | `k_m_aa * S_aa / (K_S_aa + S_aa) * X_aa * I6` | Amino acids are consumed to grow `X_aa` and form propionate, ethanol, lactate, acetate, inorganic carbon, and inorganic nitrogen. |
| `Uptake of LCFA` | `k_m_fa * S_fa / (K_S_fa + S_fa) * X_fa * I7` | LCFA are consumed to grow `X_fa` and form propionate, ethanol, lactate, acetate, carbon, and nitrogen demand. |
| `Uptake of acetate_et` | `k_m_ac * S_ac * S_et / (K_S_ac*S_ac + K_S_ac_et*S_et + S_ac*S_et + eps) * X_ac_et * I11` | Acetate and ethanol are converted by `X_ac_et`, producing butyrate, hydrogen, carbon correction, and new biomass. |
| `Uptake of acetate_lac` | `k_m_ac * S_ac * S_lac / (K_S_ac*S_ac + K_S_ac_lac*S_lac + S_ac*S_lac + eps) * X_ac_lac * I11` | Acetate and lactate are converted by `X_ac_lac`, producing butyrate, hydrogen, carbon correction, and new biomass. |
| `Uptake of propionate_et` | `k_m_pro * S_pro * S_et / (K_S_pro*S_pro + K_S_pro_et*S_et + S_pro*S_et + eps) * X_chain_et * I10` | Propionate and ethanol are converted by chain-elongating biomass, producing valerate, hydrogen, carbon correction, and `X_chain_et`. |
| `Uptake of propionate_lac` | `k_m_pro * S_pro * S_lac / (K_S_pro*S_pro + K_S_pro_lac*S_lac + S_pro*S_lac + eps) * X_chain_lac * I10` | Propionate and lactate are converted by chain-elongating biomass, producing valerate, hydrogen, carbon correction, and `X_chain_lac`. |
| `Uptake of butyrate_et` | `k_m_bu * S_bu * S_et / (K_S_bu*S_bu + K_S_bu_et*S_et + S_bu*S_et + eps) * X_chain_et * I14` | Butyrate and ethanol are converted to caproate, hydrogen, carbon correction, and `X_chain_et`. |
| `Uptake of butyrate_lac` | `k_m_bu * S_bu * S_lac / (K_S_bu*S_bu + K_S_bu_lac*S_lac + S_bu*S_lac + eps) * X_chain_lac * I14` | Butyrate and lactate are converted to caproate, hydrogen, carbon correction, and `X_chain_lac`. |
| `Uptake of valerate` | `k_m_va * S_va / (K_S_va + S_va) * X_VFA_deg * I15` | Valerate is degraded to propionate and `X_VFA_deg`. |
| `Uptake of caproate` | `k_m_cap * S_cap / (K_S_cap + S_cap) * X_VFA_deg * I13` | Caproate is degraded to acetate and `X_VFA_deg`. |
| `Uptake of butyrate` | `k_m_bu_deg * S_bu / (K_S_bu + S_bu) * X_VFA_deg * I13` | Butyrate is degraded to acetate and `X_VFA_deg`. |
| `Methanogenessis from acetate and h2` | `k_m_h2_Me_ac * S_gas_h2 * S_ac / (K_S_h2_Me_ac*S_gas_h2 + K_S_ac_Me*S_ac + S_ac*S_gas_h2 + eps) * X_Me_ac * I12` | Acetate and gas-phase hydrogen are converted to dissolved methane `S_ch4`, methanogenic biomass `X_Me_ac`, carbon correction, and nitrogen demand. |
| `Methanogenessis from CO2 and h2` | `k_m_h2_Me_CO2 * S_gas_h2 * S_gas_co2 / (K_S_h2_Me_CO2*S_gas_h2 + K_S_CO2_Me*S_gas_co2 + S_gas_co2*S_gas_h2 + eps) * X_Me_CO2 * I12` | Gas-phase hydrogen and carbon dioxide are converted to gas-phase methane `S_gas_ch4`, `X_Me_CO2`, gas-phase carbon dioxide change, and nitrogen demand. |
| `Uptake of ethanol` | `k_m_et * S_et / (K_S_et + S_et) * X_et * I16` | Ethanol is oxidized to acetate, inorganic carbon, and `X_et`. |
| `Uptake of lactate` | `k_m_lac * S_lac / (K_S_lac + S_lac) * X_lac * I16` | Lactate is oxidized to propionate, inorganic carbon, and `X_lac`. |
| `Decay of Xsu` | `k_dec_X_su * X_su` | Sugar degraders decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of Xaa` | `k_dec_X_aa * X_aa` | Amino acid degraders decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of Xfa` | `k_dec_X_fa * X_fa` | LCFA degraders decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_ac_et` | `k_dec_X_ac * X_ac_et` | Acetate-ethanol consumers decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_ac_lac` | `k_dec_X_ac * X_ac_lac` | Acetate-lactate consumers decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of Xpro` | Not currently assigned in `e_adm_ode_sys` | Listed in the reaction set but not given a rate in the current ODE block. If kept, it behaves as zero-flux unless populated elsewhere. |
| `Decay of X_chain_et` | `k_dec_X_chain_et * X_chain_et` | Ethanol-linked chain elongators decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_chain_lac` | `k_dec_X_chain_lac * X_chain_lac` | Lactate-linked chain elongators decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_VFA_deg` | `k_dec_X_VFA_deg * X_VFA_deg` | VFA degraders decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_Me_ac` | `k_dec_X_Me_ac * X_Me_ac` | Acetate/hydrogen methanogens decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of X_Me_CO2` | `k_dec_X_Me_CO2 * X_Me_CO2` | CO2/hydrogen methanogens decay to TSS, returning biomass carbon and nitrogen. |
| `Decay of Xet` | `k_dec_X_et * X_et` | Ethanol oxidizers decay. This rate is present in the ODE; the current stoichiometric builder should be checked if this flux is expected to affect states. |
| `Decay of Xlac` | `k_dec_X_lac * X_lac` | Lactate oxidizers decay. This rate is present in the ODE; the current stoichiometric builder should be checked if this flux is expected to affect states. |
| `Acid Base Equilibrium (Va)` | `k_A_B_va * (S_va_ion * (K_a_va + S_H) - K_a_va * S_va)` | Valerate speciation between `S_va` and `S_va_ion`. |
| `Acid Base Equilibrium (Bu)` | `k_A_B_bu * (S_bu_ion * (K_a_bu + S_H) - K_a_bu * S_bu)` | Butyrate speciation between `S_bu` and `S_bu_ion`. |
| `Acid Base Equilibrium (Pro)` | `k_A_B_pro * (S_pro_ion * (K_a_pro + S_H) - K_a_pro * S_pro)` | Propionate speciation between `S_pro` and `S_pro_ion`. |
| `Acid Base Equilibrium (Cap)` | `k_A_B_cap * (S_cap_ion * (K_a_cap + S_H) - K_a_cap * S_cap)` | Caproate speciation between `S_cap` and `S_cap_ion`. |
| `Acid Base Equilibrium (Lac)` | `k_A_B_lac * (S_lac_ion * (K_a_lac + S_H) - K_a_lac * S_lac)` | Lactate speciation between `S_lac` and `S_lac_ion`. |
| `Acid Base Equilibrium (Ac)` | `k_A_B_ac * (S_ac_ion * (K_a_ac + S_H) - K_a_ac * S_ac)` | Acetate speciation between `S_ac` and `S_ac_ion`. |
| `Acid Base Equilibrium (CO2)` | `k_A_B_co2 * (S_hco3_ion * (K_a_co2 + S_H) - K_a_co2 * S_IC)` | Carbonate speciation through bicarbonate `S_hco3_ion` and inorganic carbon `S_IC`. |
| `Acid Base Equilibrium (In)` | `k_A_B_IN * (S_nh3 * (K_a_IN + S_H) - K_a_IN * S_IN)` | Nitrogen speciation between ammonia `S_nh3`, ammonium `S_nh4_ion`, and total inorganic nitrogen `S_IN`. |
| `Gas Transfer H2` | `max(0, k_L_a * (S_h2 - 16 * K_H_h2 * p_gas_h2))` | Dissolved hydrogen leaves the liquid and enters `S_gas_h2` with `V_liq / V_gas` scaling. |
| `Gas Transfer CH4` | `max(0, k_L_a * (S_ch4 - 64 * K_H_ch4 * p_gas_ch4))` | Dissolved methane leaves the liquid and enters `S_gas_ch4` with `V_liq / V_gas` scaling. |
| `Gas Transfer CO2` | `max(0, k_L_a * (S_co2 - K_H_co2 * p_gas_co2))` | Dissolved carbon dioxide leaves the liquid and enters `S_gas_co2` with `V_liq / V_gas` scaling. |

### e-ADM Inhibition Factors

| Factor | Formula pattern | Applied to |
| --- | --- | --- |
| `I5` | `I_pH_aa * I_IN_lim` | Sugar uptake. |
| `I6` | `I5` | Amino acid uptake. |
| `I7` | `I_pH_aa * I_IN_lim * I_h2_fa` | LCFA uptake. |
| `I10` | `I_pH_pro * I_IN_lim * I_h2_pro` | Propionate plus ethanol or lactate uptake. |
| `I11` | `I_pH_ac * I_IN_lim * I_nh3` | Acetate plus ethanol or lactate uptake. |
| `I12` | `I_pH_h2 * I_IN_lim` | Methanogenesis. |
| `I13` | `I_pH_cap * I_IN_lim * I_h2_c4` | Caproate uptake and direct butyrate degradation. |
| `I14` | `I_pH_bu * I_IN_lim * I_h2_c4` | Butyrate plus ethanol or lactate uptake. |
| `I15` | `I_pH_va * I_IN_lim * I_h2_c4` | Valerate uptake. |
| `I16` | `I_IN_lim * I_nh3 * I_pH_aa * I_h2_oxidation` | Ethanol and lactate oxidation. |

## Methane-Producing Reactions

| Model | Methane-producing reactions |
| --- | --- |
| ADM1 | `Uptake of acetate`, `Uptake of Hydrogen` |
| e-ADM | `Methanogenessis from acetate and h2`, `Methanogenessis from CO2 and h2` |

`Gas Transfer CH4` does not create methane. It moves methane between dissolved `S_ch4` and gas-phase `S_gas_ch4`.

## State Vectors

### ADM1 Species

`S_su`, `S_aa`, `S_fa`, `S_va`, `S_bu`, `S_pro`, `S_ac`, `S_h2`, `S_ch4`, `S_IC`, `S_IN`, `S_I`, `X_xc`, `X_ch`, `X_pr`, `X_li`, `X_su`, `X_aa`, `X_fa`, `X_c4`, `X_pro`, `X_ac`, `X_h2`, `X_I`, `S_cation`, `S_anion`, `S_H_ion`, `S_va_ion`, `S_bu_ion`, `S_pro_ion`, `S_ac_ion`, `S_hco3_ion`, `S_co2`, `S_nh3`, `S_nh4_ion`, `S_gas_h2`, `S_gas_ch4`, `S_gas_co2`

### e-ADM Species

`S_su`, `S_aa`, `S_fa`, `S_va`, `S_bu`, `S_pro`, `S_et`, `S_lac`, `S_cap`, `S_ac`, `S_h2`, `S_ch4`, `S_IC`, `S_IN`, `S_I`, `TSS`, `TDS`, `X_ch`, `X_pr`, `X_li`, `X_su`, `X_aa`, `X_ac_et`, `X_ac_lac`, `X_fa`, `X_VFA_deg`, `X_et`, `X_lac`, `X_chain_et`, `X_chain_lac`, `X_Me_ac`, `X_Me_CO2`, `X_I`, `S_cation`, `S_anion`, `S_H_ion`, `S_va_ion`, `S_bu_ion`, `S_pro_ion`, `S_cap_ion`, `S_lac_ion`, `S_ac_ion`, `S_hco3_ion`, `S_co2`, `S_nh3`, `S_nh4_ion`, `S_gas_h2`, `S_gas_ch4`, `S_gas_co2`

## Minimal CLI Run

```bash
adtoolbox ADM adm1 --models-json /path/to/ADToolbox/models.json --report csv
```

```bash
adtoolbox ADM e-adm --models-json /path/to/ADToolbox/models.json --report csv
```
