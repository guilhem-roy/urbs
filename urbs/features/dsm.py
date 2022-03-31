import math
import pyomo.core as pyomo
import pandas as pd
import xarray as xr


def add_dsm(m, m_parameters):

    # modelled Demand Side Management time steps (downshift):
    # downshift effective in tt to compensate for upshift in t
    m_parameters['tt'] = pd.Index(m_parameters['timesteps'][1:], name="tt")
    # m.tt = pyomo.Set(
    #     within=m.t,
    #     initialize=m.timesteps[1:],
    #     ordered=True,
    #     doc='Set of additional DSM time steps')

    # DSM Tuples
    m_parameters['dsm_site_tuples'] = pd.Index(m_parameters['dsm_dict']["delay"].keys(),
            name="dsm_site_tuples", tupleize_cols=False)
    # m.dsm_site_tuples = pyomo.Set(
    #     within=m.stf*m.sit*m.com,
    #     initialize=tuple(m.dsm_dict["delay"].keys()),
    #     doc='Combinations of possible dsm by site, e.g. '
    #         '(2020, Mid, Elec)')
    m_parameters['dsm_down_tuples'] = pd.Index([(t, tt, stf, site, commodity)
        for (t, tt, stf, site, commodity)
        in dsm_down_time_tuples(m_parameters['timesteps'][1:],
            m_parameters['dsm_site_tuples'],
            m_parameters)],
        name="dsm_down_tuples", tupleize_cols=False)
    # m.dsm_down_tuples = pyomo.Set(
    #     within=m.tm*m.tm*m.stf*m.sit*m.com,
    #     initialize=[(t, tt, stf, site, commodity)
    #                 for (t, tt, stf, site, commodity)
    #                 in dsm_down_time_tuples(m.timesteps[1:],
    #                                         m.dsm_site_tuples,
    #                                         m)],
    #     doc='Combinations of possible dsm_down combinations, e.g. '
    #         '(5001,5003,2020,Mid,Elec)')

    # Variables
    m.add_variables(name="dsm_up", coords=[m_parameters['tm'],
        m_parameters['dsm_site_tuples']], lower=0)
    # m.dsm_up = pyomo.Var(
    #     m.tm, m.dsm_site_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='DSM upshift')
    m.add_variables(name="dsm_down", coords=[m_parameters['dsm_down_tuples']],
            lower=0)
    # m.dsm_down = pyomo.Var(
    #     m.dsm_down_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='DSM downshift')

    # TODO : Improve preformance, avoid loops
    # XXX : Remove this return
    return
    # DSM rules
    for (stf, sit, com) in m_parameters['dsm_site_tuples']:
        delay = max(int(1 / m_parameters['dt'] *
            m_parameters['dsm_dict']['delay'][(stf, sit, com)]), 1)
        for tm in m_parameters['tm']:
            dsm_down_sum = m.variables["dsm_down"].loc[[(tm, tt, stf, sit, com)
                for tt in dsm_time_tuples(tm, m_parameters['timesteps'][1:], delay
                    )]].sum()
            m.add_constraints(dsm_down_sum -
                    m.variables['dsm_up'].loc[tm, (stf, sit, com)].sum() *
                    m_parameters['dsm_dict']['eff'][(stf, sit, com)], "=", 0,
                    name="def_dsm_variables"+str((tm, stf, sit, com)))
            # m.def_dsm_variables = pyomo.Constraint(
            #     m.tm, m.dsm_site_tuples,
            #     rule=def_dsm_variables_rule,
            #     doc='DSMup * efficiency factor n == DSMdo (summed)')

    m.add_constraints(m.variables['dsm_up'], "<=", 
            m_parameters['dt'] * 
            xr.DataArray.from_series(pd.Series(
                m_parameters['dsm_dict']['cap-max-up'],
                name="dsm_site_tuples")),
            name="res_dsm_upward")
    # m.res_dsm_upward = pyomo.Constraint(
    #     m.tm, m.dsm_site_tuples,
    #     rule=res_dsm_upward_rule,
    #     doc='DSMup <= Cup (threshold capacity of DSMup)')

    for (stf, sit, com) in m_parameters['dsm_site_tuples']:
        delay =  max(int(1 / m_parameters['dt'] *
            m_parameters['dsm_dict']['delay'][(stf, sit, com)]), 1)
        for tm in m_parameters['tm']:
            dsm_down_sum = m.variables["dsm_down"].loc[[(t, tm, stf, sit, com)
                for t in dsm_time_tuples(tm, m_parameters['timesteps'][1:],
                    delay)]].sum()
            m.add_constraints(dsm_down_sum, "<=", m_parameters['dt'] *
                    m_parameters['dsm_dict']['cap-max-do'][(stf, sit, com)],
                    name="res_dsm_downward"+str((tm, stf, sit, com)))
            # m.res_dsm_downward = pyomo.Constraint(
            #     m.tm, m.dsm_site_tuples,
            #     rule=res_dsm_downward_rule,
            #     doc='DSMdo (summed) <= Cdo (threshold capacity of DSMdo)')
            max_dsm_limit = (m_parameters['dt'] *
                    max(m_parameters['dsm_dict']['cap-max-up'][(stf, sit, com)],
                        m_parameters['dsm_dict']['cap-max-do'][(stf, sit, com)]))
            m.add_constraints(m.variables['dsm_up'].loc[tm, (stf, sit, com)].sum() +
                    dsm_down_sum, "<=", max_dsm_limit,
                    name="res_dsm_maximum"+str((tm, stf, sit, com)))
            # m.res_dsm_maximum = pyomo.Constraint(
            #     m.tm, m.dsm_site_tuples,
            #     rule=res_dsm_maximum_rule,
            #     doc='DSMup + DSMdo (summed) <= max(Cup,Cdo)')

    for (stf, sit, com) in m_parameters['dsm_site_tuples']:
        recov = max(int(1 / m_parameters['dt'] *
            m_parameters['dsm_dict']['recov'][(stf, sit, com)]), 1)
        for tm in m_parameters['tm']:
            dsm_up_sum = m.variables['dsm_up'].loc[dsm_recovery(tm,
                m_parameters['timesteps'][1:], recov), (stf, sit, com)].sum()
            m.add_constraints(dsm_up_sum, "<=",
                    m_parameters['dsm_dict']['cap-max-up'][(stf, sit, com)] *
                    m_parameters['dsm_dict']['delay'][(stf, sit, com)],
                    name="res_dsm_discovery"+str((tm, stf, sit, com)))
            # m.res_dsm_recovery = pyomo.Constraint(
            #     m.tm, m.dsm_site_tuples,
            #     rule=res_dsm_recovery_rule,
            #     doc='DSMup(t, t + recovery time R) <= Cup * delay time L')

    return m


# demand side management (DSM) constraints

# DSMup == DSMdo * efficiency factor n
def def_dsm_variables_rule(m, tm, stf, sit, com):
    dsm_down_sum = 0
    for tt in dsm_time_tuples(tm,
                              m.timesteps[1:],
                              max(int(1 / m.dt *
                                  m.dsm_dict['delay'][(stf, sit, com)]), 1)):
        dsm_down_sum += m.dsm_down[tm, tt, stf, sit, com]
    return dsm_down_sum == (m.dsm_up[tm, stf, sit, com] *
                            m.dsm_dict['eff'][(stf, sit, com)])


# DSMup <= Cup (threshold capacity of DSMup)
def res_dsm_upward_rule(m, tm, stf, sit, com):
    return m.dsm_up[tm, stf, sit, com] <= (m.dt *
                                           m.dsm_dict['cap-max-up']
                                           [(stf, sit, com)])


# DSMdo <= Cdo (threshold capacity of DSMdo)
def res_dsm_downward_rule(m, tm, stf, sit, com):
    dsm_down_sum = 0
    for t in dsm_time_tuples(tm,
                             m.timesteps[1:],
                             max(int(1 / m.dt *
                                 m.dsm_dict['delay'][(stf, sit, com)]), 1)):
        dsm_down_sum += m.dsm_down[t, tm, stf, sit, com]
    return dsm_down_sum <= (m.dt * m.dsm_dict['cap-max-do'][(stf, sit, com)])


# DSMup + DSMdo <= max(Cup,Cdo)
def res_dsm_maximum_rule(m, tm, stf, sit, com):
    dsm_down_sum = 0
    for t in dsm_time_tuples(tm,
                             m.timesteps[1:],
                             max(int(1 / m.dt *
                                 m.dsm_dict['delay'][(stf, sit, com)]), 1)):
        dsm_down_sum += m.dsm_down[t, tm, stf, sit, com]

    max_dsm_limit = m.dt * max(m.dsm_dict['cap-max-up'][(stf, sit, com)],
                               m.dsm_dict['cap-max-do'][(stf, sit, com)])
    return m.dsm_up[tm, stf, sit, com] + dsm_down_sum <= max_dsm_limit


# DSMup(t, t + recovery time R) <= Cup * delay time L
def res_dsm_recovery_rule(m, tm, stf, sit, com):
    dsm_up_sum = 0
    for t in dsm_recovery(tm,
                          m.timesteps[1:],
                          max(int(1 / m.dt *
                              m.dsm_dict['recov'][(stf, sit, com)]), 1)):
        dsm_up_sum += m.dsm_up[t, stf, sit, com]
    return dsm_up_sum <= (m.dsm_dict['cap-max-up'][(stf, sit, com)] *
                          m.dsm_dict['delay'][(stf, sit, com)])


# DSM surplus
def dsm_surplus(m, tm, stf, sit, com):
    """ called in vertex rule
        calculate dsm surplus"""
    if (stf, sit, com) in m.dsm_site_tuples:
        return (- m.dsm_up[tm, stf, sit, com] +
                sum(m.dsm_down[t, tm, stf, sit, com]
                    for t in dsm_time_tuples(
                    tm, m.timesteps[1:],
                    max(int(1 / m.dt *
                        m.dsm_dict['delay'][(stf, sit, com)]), 1))))
    else:
        return 0


def dsm_down_time_tuples(time, sit_com_tuple, m_parameters):
    """ Dictionary for the two time instances of DSM_down
    Args:
        time: list with time indices
        sit_com_tuple: a list of (site, commodity) tuples
        m_parameters: model parameters
    Returns:
        A list of possible time tuples depending on site and commodity
    """

    delay = m_parameters['dsm_dict']['delay']
    ub = max(time)
    lb = min(time)
    time_list = []

    for (stf, site, commodity) in sit_com_tuple:
        for step1 in time:
            for step2 in range(step1 -
                               max(int(delay[stf, site, commodity] /
                                   m_parameters['dt']), 1),
                               step1 +
                               max(int(delay[stf, site, commodity] /
                                   m_parameters['dt']), 1) + 1):
                if lb <= step2 <= ub:
                    time_list.append((step1, step2, stf, site, commodity))

    return time_list


def dsm_time_tuples(timestep, time, delay):
    """ Tuples for the two time instances of DSM_down
    Args:
        timestep: current timestep
        time: list with time indices
        delay: allowed dsm delay in particular site and commodity
    Returns:
        A list of possible time tuples of a current time step in a specific
        stf, site and commodity
    """

    ub = max(time)
    lb = min(time)

    time_list = list()

    for step in range(timestep - delay, timestep + delay + 1):
        if step >= lb and step <= ub:
            time_list.append(step)

    return time_list


def dsm_recovery(timestep, time, recov):
    """ Time frame for the allowed time indices in case of recovery
    Args:
        timestep: current timestep
        time: list with time indices
        recov: allowed dsm recovery in particular site and commodity
    Returns:
        A list of possible time indices which are within the modelled time
        area
    """

    ub = max(time)

    time_list = list()

    for step in range(timestep, timestep+recov):
        if step <= ub:
            time_list.append(step)

    return time_list
