from pybamm import exp, tanh, constants, Parameter, FunctionParameter, ParameterValues


def electrolyte_conductivity_Ecker2015(c_e, T):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration [1, 2, 3].

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    # mol/m^3 to mol/l
    cm = 1e-3 * c_e

    # value at T = 296K
    sigma_e_296 = 0.2667 * cm**3 - 1.2983 * cm**2 + 1.7919 * cm + 0.1726

    # add temperature dependence
    E_k_e = 1.71e4
    C = 296 * exp(E_k_e / (constants.R * 296))
    sigma_e = C * sigma_e_296 * exp(-E_k_e / (constants.R * T)) / T

    return sigma_e


def electrolyte_diffusivity_Ecker2015(c_e, T):
    """
    Diffusivity of LiPF6 in EC:DMC as a function of ion concentration [1, 2, 3].

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    # The diffusivity epends on the electrolyte conductivity
    inputs = {"Electrolyte concentration [mol.m-3]": c_e, "Temperature [K]": T}
    sigma_e = FunctionParameter("Electrolyte conductivity [S.m-1]", inputs)

    D_c_e = (constants.k_b / (constants.F * constants.q_e)) * sigma_e * T / c_e

    return D_c_e


def graphite_diffusivity_Ecker2015(sto, T):
    """
    Graphite diffusivity as a function of stochiometry [1, 2, 3].

    References
    ----------
     .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_ref = 8.4e-13 * exp(-11.3 * sto) + 8.2e-15
    E_D_s = 3.03e4
    arrhenius = exp(-E_D_s / (constants.R * T)) * exp(E_D_s / (constants.R * 296))

    return D_ref * arrhenius


def graphite_electrolyte_exchange_current_density_Ecker2015(c_e, c_s_surf, T):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """

    k_ref = 1.995 * 1e-10

    # multiply by Faraday's constant to get correct units
    m_ref = constants.F * k_ref  # (A/m2)(mol/m3)**1.5 - includes ref concentrations
    E_r = 53400

    arrhenius = exp(-E_r / (constants.R * T)) * exp(E_r / (constants.R * 296.15))

    c_n_max = Parameter("Maximum concentration in negative electrode [mol.m-3]")

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_n_max - c_s_surf) ** 0.5
    )


def graphite_ocp_Ecker2015_function(sto):
    """
    Graphite OCP as a function of stochiometry [1, 2, 3].

    References
    ----------
     .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """

    # Graphite negative electrode from Ecker, Kabitz, Laresgoiti et al.
    # Analytical fit (WebPlotDigitizer + gnuplot)
    a = 0.716502
    b = 369.028
    c = 0.12193
    d = 35.6478
    e = 0.0530947
    g = 0.0169644
    h = 27.1365
    i = 0.312832
    j = 0.0199313
    k = 28.5697
    m = 0.614221
    n = 0.931153
    o = 36.328
    p = 1.10743
    q = 0.140031
    r = 0.0189193
    s = 21.1967
    t = 0.196176

    u_eq = (
        a * exp(-b * sto)
        + c * exp(-d * (sto - e))
        - r * tanh(s * (sto - t))
        - g * tanh(h * (sto - i))
        - j * tanh(k * (sto - m))
        - n * exp(o * (sto - p))
        + q
    )

    return u_eq


def nco_diffusivity_Ecker2015(sto, T):
    """
    NCO diffusivity as a function of stochiometry [1, 2, 3].

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_ref = 3.7e-13 - 3.4e-13 * exp(-12 * (sto - 0.62) * (sto - 0.62))
    E_D_s = 8.06e4
    arrhenius = exp(-E_D_s / (constants.R * T)) * exp(E_D_s / (constants.R * 296.15))

    return D_ref * arrhenius


def nco_electrolyte_exchange_current_density_Ecker2015(c_e, c_s_surf, T):
    """
    Exchange-current density for Butler-Volmer reactions between NCO and LiPF6 in
    EC:DMC [1, 2, 3].

    References
    ----------
       .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """

    k_ref = 5.196e-11

    # multiply by Faraday's constant to get correct units
    m_ref = constants.F * k_ref  # (A/m2)(mol/m3)**1.5 - includes ref concentrations

    E_r = 4.36e4
    arrhenius = exp(-E_r / (constants.R * T)) * exp(E_r / (constants.R * 296.15))

    c_p_max = Parameter("Maximum concentration in positive electrode [mol.m-3]")

    return (
        m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_p_max - c_s_surf) ** 0.5
    )


def nco_ocp_Ecker2015_function(sto):
    """
    NCO OCP as a function of stochiometry [1, 2, 3].

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
       Stochiometry of material (li-fraction)

    """

    # LiNiCo from Ecker, Kabitz, Laresgoiti et al.
    # Analytical fit (WebPlotDigitizer + gnuplot)
    a = -2.35211
    c = 0.0747061
    d = 31.886
    e = 0.0219921
    g = 0.640243
    h = 5.48623
    i = 0.439245
    j = 3.82383
    k = 4.12167
    m = 0.176187
    n = 0.0542123
    o = 18.2919
    p = 0.762272
    q = 4.23285
    r = -6.34984
    s = 2.66395
    t = 0.174352

    u_eq = (
        a * sto
        - c * tanh(d * (sto - e))
        - r * tanh(s * (sto - t))
        - g * tanh(h * (sto - i))
        - j * tanh(k * (sto - m))
        - n * tanh(o * (sto - p))
        + q
    )
    return u_eq


parameter_values = ParameterValues(
    {
        "1 + dlnf/dlnc": 1.0,
        "Ambient temperature [K]": 298.15,
        "Cation transference number": 0.26,
        "Cell cooling surface area [m2]": 0.0172,
        "Cell volume [m3]": 1.52e-06,
        "Current function [A]": 0.15652,
        "Electrode height [m]": 0.101,
        "Electrode width [m]": 0.085,
        "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Ecker2015,
        "Electrolyte diffusivity [m2.s-1]": electrolyte_diffusivity_Ecker2015,
        "Initial concentration in electrolyte [mol.m-3]": 1000.0,
        "Initial concentration in negative electrode [mol.m-3]": 26120.05,
        "Initial concentration in positive electrode [mol.m-3]": 12630.8,
        "Initial temperature [K]": 298.15,
        "Lower voltage cut-off [V]": 2.5,
        "Maximum concentration in negative electrode [mol.m-3]": 31920.0,
        "Maximum concentration in positive electrode [mol.m-3]": 48580.0,
        "Negative current collector conductivity [S.m-1]": 58411000.0,
        "Negative current collector density [kg.m-3]": 8933.0,
        "Negative current collector specific heat capacity [J.kg-1.K-1]": 385.0,
        "Negative current collector thermal conductivity [W.m-1.K-1]": 398.0,
        "Negative current collector thickness [m]": 1.4e-05,
        "Negative electrode Bruggeman coefficient (electrode)": 0.0,
        "Negative electrode Bruggeman coefficient (electrolyte)": 1.6372789338386007,
        "Negative electrode OCP [V]": graphite_ocp_Ecker2015_function,
        "Negative electrode OCP entropic change [V.K-1]": 0.0,
        "Negative electrode active material volume fraction": 0.372403,
        "Negative electrode conductivity [S.m-1]": 14.0,
        "Negative electrode density [kg.m-3]": 1555.0,
        "Negative electrode diffusivity [m2.s-1]": graphite_diffusivity_Ecker2015,
        "Negative electrode electrons in reaction": 1.0,
        "Negative electrode exchange-current density [A.m-2]": graphite_electrolyte_exchange_current_density_Ecker2015,
        "Negative electrode porosity": 0.329,
        "Negative electrode specific heat capacity [J.kg-1.K-1]": 1437.0,
        "Negative electrode thermal conductivity [W.m-1.K-1]": 1.58,
        "Negative electrode thickness [m]": 7.4e-05,
        "Negative particle radius [m]": 1.37e-05,
        "Nominal cell capacity [A.h]": 0.15625,
        "Number of cells connected in series to make a battery": 1.0,
        "Number of electrodes connected in parallel to make a cell": 1.0,
        "Positive current collector conductivity [S.m-1]": 36914000.0,
        "Positive current collector density [kg.m-3]": 2702.0,
        "Positive current collector specific heat capacity [J.kg-1.K-1]": 903.0,
        "Positive current collector thermal conductivity [W.m-1.K-1]": 238.0,
        "Positive current collector thickness [m]": 1.5e-05,
        "Positive electrode Bruggeman coefficient (electrode)": 0.0,
        "Positive electrode Bruggeman coefficient (electrolyte)": 1.5442267190786427,
        "Positive electrode OCP [V]": nco_ocp_Ecker2015_function,
        "Positive electrode OCP entropic change [V.K-1]": 0.0,
        "Positive electrode active material volume fraction": 0.40832,
        "Positive electrode conductivity [S.m-1]": 68.1,
        "Positive electrode density [kg.m-3]": 2895.0,
        "Positive electrode diffusivity [m2.s-1]": nco_diffusivity_Ecker2015,
        "Positive electrode electrons in reaction": 1.0,
        "Positive electrode exchange-current density [A.m-2]": nco_electrolyte_exchange_current_density_Ecker2015,
        "Positive electrode porosity": 0.296,
        "Positive electrode specific heat capacity [J.kg-1.K-1]": 1270.0,
        "Positive electrode thermal conductivity [W.m-1.K-1]": 1.04,
        "Positive electrode thickness [m]": 5.4e-05,
        "Positive particle radius [m]": 6.5e-06,
        "Reference temperature [K]": 296.15,
        "Separator Bruggeman coefficient (electrolyte)": 1.9804586773134945,
        "Separator density [kg.m-3]": 1017.0,
        "Separator porosity": 0.508,
        "Separator specific heat capacity [J.kg-1.K-1]": 1978.0,
        "Separator thermal conductivity [W.m-1.K-1]": 0.34,
        "Separator thickness [m]": 2e-05,
        "Total heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Typical current [A]": 0.15652,
        "Typical electrolyte concentration [mol.m-3]": 1000.0,
        "Upper voltage cut-off [V]": 4.2,
    }
)
