import json
import time
from pathlib import Path

import numpy as np
import pytest

from adtoolbox import adm, utils


DATABASE_ROOT = Path("/Users/parsaghadermarzi/Desktop/Academics/Projects/Database/ADToolbox")


def _load_parameter_set(folder, prefix):
    base = DATABASE_ROOT / folder
    if not base.exists():
        pytest.skip(f"ADToolbox database not available at {DATABASE_ROOT}")
    return {
        suffix: json.loads((base / f"{prefix}_{suffix}.json").read_text())
        for suffix in [
            "species",
            "reactions",
            "model_parameters",
            "base_parameters",
            "initial_conditions",
            "inlet_conditions",
        ]
    }


def _model_from_params(params, ode_system, build_stoichiometric_matrix):
    return adm.Model(
        model_parameters=params["model_parameters"],
        base_parameters=params["base_parameters"],
        initial_conditions=params["initial_conditions"],
        inlet_conditions=params["inlet_conditions"],
        feed=adm.DEFAULT_FEED,
        reactions=params["reactions"],
        species=params["species"],
        ode_system=ode_system,
        build_stoichiometric_matrix=build_stoichiometric_matrix,
    )


def test_adm1_methane_is_produced_by_acetate_and_hydrogen_uptake():
    params = _load_parameter_set("adm1", "adm1")
    S = adm.build_adm1_stoichiometric_matrix(
        params["base_parameters"],
        params["model_parameters"],
        params["reactions"],
        params["species"],
        adm.DEFAULT_FEED,
    )

    methane_row = params["species"].index("S_ch4")
    producers = {
        params["reactions"][idx]
        for idx in np.where(S[methane_row, :] > 0)[0]
    }

    assert producers == {"Uptake of acetate", "Uptake of Hydrogen"}


def test_e_adm_methanogenesis_produces_methane_and_gas_transfer_uses_adm1_scaling():
    params = _load_parameter_set("e_adm_2", "e_adm_2")
    S = adm.build_e_adm_stoichiometric_matrix(
        params["base_parameters"],
        params["model_parameters"],
        params["reactions"],
        params["species"],
        adm.DEFAULT_FEED,
    )

    acetate_methanogenesis = params["reactions"].index("Methanogenessis from acetate and h2")
    co2_methanogenesis = params["reactions"].index("Methanogenessis from CO2 and h2")
    # Both methanogenesis routes now produce DISSOLVED S_ch4, which reaches the
    # headspace through Gas Transfer CH4 (hydrogenotrophic previously emitted
    # gas-phase CH4 directly, bypassing the liquid phase).
    assert S[params["species"].index("S_ch4"), acetate_methanogenesis] > 0
    assert S[params["species"].index("S_ch4"), co2_methanogenesis] > 0

    ch4_transfer = params["reactions"].index("Gas Transfer CH4")
    volume_ratio = params["base_parameters"]["V_liq"] / params["base_parameters"]["V_gas"]
    assert S[params["species"].index("S_ch4"), ch4_transfer] == -1
    assert S[params["species"].index("S_gas_ch4"), ch4_transfer] == volume_ratio


def test_monod_limitation_increases_with_substrate():
    half_saturation = 0.01
    low = adm._monod_limitation(0.001, half_saturation)
    high = adm._monod_limitation(1.0, half_saturation)

    assert 0 < low < high < 1


def test_all_models_json_loads_model_by_key(tmp_path):
    params = _load_parameter_set("adm1", "adm1")
    models_json = tmp_path / "models.json"
    models_json.write_text(json.dumps({"adm1": params, "e_adm": {"species": []}}))

    loaded = utils.load_model_json(str(models_json), "adm1")

    assert loaded["species"] == params["species"]
    assert loaded["reactions"] == params["reactions"]


def test_e_adm_nitrogen_acid_base_rate_uses_s_in_not_s_ic():
    params = _load_parameter_set("e_adm_2", "e_adm_2")
    model = _model_from_params(
        params,
        adm.e_adm_ode_sys,
        adm.build_e_adm_stoichiometric_matrix,
    )
    model.info = {"Fluxes": []}
    model._be_time = time.time()

    c = model.initial_conditions[:, 0].copy()
    c[model.species.index("S_IC")] = 10.0
    c[model.species.index("S_IN")] = 0.2
    c[model.species.index("S_H_ion")] = 1e-7
    h_at_flux_calculation = c[model.species.index("S_H_ion")]

    adm.e_adm_ode_sys(0, c, model)
    flux = model.info["Fluxes"][model.reactions.index("Acid Base Equilibrium (In)"), 0]
    expected = model.model_parameters["k_A_B_IN"] * (
        c[model.species.index("S_nh3")]
        * (model.model_parameters["K_a_IN"] + h_at_flux_calculation)
        - model.model_parameters["K_a_IN"] * c[model.species.index("S_IN")]
    )

    assert flux == pytest.approx(expected)


def test_e_adm_reactions_conserve_cod():
    """Every biological conversion must conserve COD (sum of COD-species
    coefficients == 0). Gas-transfer reactions are excluded because gas
    concentrations are per gas-volume (V_liq/V_gas scaling), and the DAE-handled
    acid/base equilibria are excluded. Guards the methanogenesis and chain-
    elongation stoichiometry fixes against regression."""
    params = _load_parameter_set("e_adm_2", "e_adm_2")
    sp = params["species"]
    rx = params["reactions"]
    S = adm.build_e_adm_stoichiometric_matrix(
        params["base_parameters"], params["model_parameters"], rx, sp, adm.DEFAULT_FEED,
    )
    non_cod = {"S_IC", "S_IN", "S_cation", "S_anion", "S_H_ion", "S_hco3_ion",
               "S_co2", "S_nh3", "S_nh4_ion", "S_gas_co2"}
    cod = np.array([0.0 if s in non_cod else 1.0 for s in sp])
    skip = {r for r in rx if r.startswith("Gas Transfer") or r.startswith("Acid Base")}
    offenders = {}
    for j, r in enumerate(rx):
        if r in skip:
            continue
        total = float((np.asarray(S)[:, j] * cod).sum())
        if abs(total) > 1e-6:
            offenders[r] = round(total, 5)
    assert not offenders, f"COD not conserved in: {offenders}"


def test_model_save_load_round_trip(tmp_path):
    """Model.save / Model.load reproduce the model exactly, including edited
    initial conditions and the resolved callables."""
    repo = Path(__file__).resolve().parent.parent
    p = utils.load_model_json(str(repo / "reference_data" / "models.json"), "e_adm")
    model = adm.Model(
        model_parameters=p["model_parameters"], base_parameters=p["base_parameters"],
        initial_conditions=p["initial_conditions"], inlet_conditions=p["inlet_conditions"],
        feed=adm.DEFAULT_FEED, reactions=p["reactions"], species=p["species"],
        ode_system=adm.e_adm_ode_sys,
        build_stoichiometric_matrix=adm.build_e_adm_stoichiometric_matrix,
        control_state={"S_H_ion": 10 ** -6.5}, name="e-ADM", switch="DAE",
    )
    model.update_parameters(model_parameters={"k_m_su": 42.0},
                            initial_conditions={"S_ac": 0.5})

    out = tmp_path / "model.json"
    model.save(out)
    loaded = adm.Model.load(out)

    assert loaded.name == model.name and loaded.switch == model.switch
    assert loaded.model_parameters["k_m_su"] == 42.0
    assert float(loaded.initial_conditions[loaded.species.index("S_ac"), 0]) == 0.5
    assert loaded.ode_system is adm.e_adm_ode_sys
    assert loaded.build_stoichiometric_matrix is adm.build_e_adm_stoichiometric_matrix

    t = np.linspace(0, 10, 50)
    assert np.array_equal(model.solve_model(t).y, loaded.solve_model(t).y)
