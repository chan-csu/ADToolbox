import json
import time
from pathlib import Path

import numpy as np
import pytest

from adtoolbox import adm


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
    assert S[params["species"].index("S_ch4"), acetate_methanogenesis] > 0
    assert S[params["species"].index("S_gas_ch4"), co2_methanogenesis] > 0

    ch4_transfer = params["reactions"].index("Gas Transfer CH4")
    volume_ratio = params["base_parameters"]["V_liq"] / params["base_parameters"]["V_gas"]
    assert S[params["species"].index("S_ch4"), ch4_transfer] == -1
    assert S[params["species"].index("S_gas_ch4"), ch4_transfer] == volume_ratio


def test_monod_limitation_increases_with_substrate():
    half_saturation = 0.01
    low = adm._monod_limitation(0.001, half_saturation)
    high = adm._monod_limitation(1.0, half_saturation)

    assert 0 < low < high < 1


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
