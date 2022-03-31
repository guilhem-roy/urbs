import math
import pyomo.core as pyomo
import pandas as pd
import xarray as xr


def add_storage(m, m_parameters):

    # storage (e.g. hydrogen, pump storage)
    indexlist = set()
    for key in m_parameters['storage_dict']["eff-in"]:
        indexlist.add(tuple(key)[2])
    m_parameters['sto'] = pd.Index(indexlist, name="sto")
    # m.sto = pyomo.Set(
    #     initialize=indexlist,
    #     doc='Set of storage technologies')

    # storage tuples
    m_parameters['sto_tuples'] = pd.Index(m_parameters['storage_dict']["eff-in"].keys(),
            name="sto_tuples", tupleize_cols=False)
    # m.sto_tuples = pyomo.Set(
    #     within=m.stf * m.sit * m.sto * m.com,
    #     initialize=tuple(m.storage_dict["eff-in"].keys()),
    #     doc='Combinations of possible storage by site,'
    #         'e.g. (2020,Mid,Bat,Elec)')

    # tuples for intertemporal operation
    if m_parameters['mode']['int']:
        m_parameters['operational_sto_tuples'] = pd.Index(
                [(sit, sto, com, stf, stf_later)
                for (sit, sto, com, stf, stf_later)
                in op_sto_tuples(m_parameters['sto_tuples'], m_parameters)],
                name='operational_sto_tuples', tupleize_cols=False)
        # m.operational_sto_tuples = pyomo.Set(
        #     within=m.sit * m.sto * m.com * m.stf * m.stf,
        #     initialize=[(sit, sto, com, stf, stf_later)
        #                 for (sit, sto, com, stf, stf_later)
        #                 in op_sto_tuples(m.sto_tuples, m)],
        #     doc='Processes that are still operational through stf_later'
        #         '(and the relevant years following), if built in stf'
        #         'in stf.')
        m_parameters['inst_sto_tuples'] = pd.index([(sit, sto, com, stf)
                for (sit, sto, com, stf)
                in inst_sto_tuples(m_parameters)],
                name="inst_sto_tuples", tupleize_cols=False)
        # m.inst_sto_tuples = pyomo.Set(
        #     within=m.sit * m.sto * m.com * m.stf,
        #     initialize=[(sit, sto, com, stf)
        #                 for (sit, sto, com, stf)
        #                 in inst_sto_tuples(m)],
        #     doc='Installed storages that are still operational through stf')

    # storage tuples for storages with fixed initial state
    m_parameters['sto_init_bound_tuples'] = pd.Index(
            m_parameters['stor_init_bound_dict'].keys(),
            name="sto_init_bound_tuples", tupleize_cols=False)
    # m.sto_init_bound_tuples = pyomo.Set(
    #     within=m.stf * m.sit * m.sto * m.com,
    #     initialize=tuple(m.stor_init_bound_dict.keys()),
    #     doc='storages with fixed initial state')

    # storage tuples for storages with given energy to power ratio
    m_parameters['sto_ep_ratio_tuples'] = pd.Index(m_parameters['sto_ep_ratio_dict'].keys(),
            name="sto_ep_ratio_tuples", tupleize_cols=False)
    # m.sto_ep_ratio_tuples = pyomo.Set(
    #     within=m.stf * m.sit * m.sto * m.com,
    #     initialize=tuple(m.sto_ep_ratio_dict.keys()),
    #     doc='storages with given energy to power ratio')

    # Variables
    m.add_variables(lower=0, coords=[m_parameters['sto_tuples']],
            name="cap_sto_c_new")
    # m.cap_sto_c_new = pyomo.Var(
    #     m.sto_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='New storage size (MWh)')
    m.add_variables(lower=0, coords=[m_parameters['sto_tuples']],
            name="cap_sto_p_new")
    # m.cap_sto_p_new = pyomo.Var(
    #     m.sto_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='New  storage power (MW)')

    # storage capacities as expression objects
    m_parameters['cap_sto_c'] = {
        (stf, sit, sto, com): def_storage_capacity_rule(m, m_parameters, stf, sit, sto, com)
        for (stf, sit, sto, com) in m_parameters['sto_tuples']}
    # m.cap_sto_c = pyomo.Expression(
    #     m.sto_tuples,
    #     rule=def_storage_capacity_rule,
    #     doc='Total storage size (MWh)')
    m_parameters['cap_sto_p'] = {
        (stf, sit, sto, com): def_storage_power_rule(m, m_parameters, stf, sit, sto, com)
        for (stf, sit, sto, com) in m_parameters['sto_tuples']}
    # m.cap_sto_p = pyomo.Expression(
    #     m.sto_tuples,
    #     rule=def_storage_power_rule,
    #     doc='Total storage power (MW)')
    m.add_variables(name="e_sto_in",
            coords=[m_parameters['tm'], m_parameters['sto_tuples']],
            lower=0)
    # m.e_sto_in = pyomo.Var(
    #     m.tm, m.sto_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='Power flow into storage (MW) per timestep')
    m.add_variables(name="e_sto_out",
            coords=[m_parameters['tm'], m_parameters['sto_tuples']],
            lower=0)
    # m.e_sto_out = pyomo.Var(
    #     m.tm, m.sto_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='Power flow out of storage (MW) per timestep')
    m.add_variables(name="e_sto_con",
            coords=[m_parameters['t'], m_parameters['sto_tuples']],
            lower=0)
    # Linopy doesn't like when a variable appears multiple times in
    # one constraint: it create conflicts for dimension names. Hence the need
    # to have two variables to represent e_sto_con
    # m.e_sto_con = pyomo.Var(
    #     m.t, m.sto_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='Energy content of storage (MWh) in timestep')

    # storage rules
    m.add_constraints(m.variables['e_sto_con'][1:,] -
        m.variables['e_sto_con'].shift({"t":1})[1:,] *
        xr.DataArray.from_series(pd.Series({s:
            (1 - m_parameters['storage_dict']['discharge'] [s]) **
            m_parameters['dt']
            for s in m_parameters['sto_tuples']}, name="sto_tuples")) -
        m.variables['e_sto_in'][1:,] *
        xr.DataArray.from_series(pd.Series(m_parameters['storage_dict']['eff-in'], name="sto_tuples")) +
        m.variables['e_sto_out'][1:,] *
        xr.DataArray.from_series(1/pd.Series(m_parameters['storage_dict']['eff-out'], name="sto_tuples")),
        "=", 0,
        name="def_storage_state")
    # m.def_storage_state = pyomo.Constraint(
    #     m.tm, m.sto_tuples,
    #     rule=def_storage_state_rule,
    #     doc='storage[t] = (1 - sd) * storage[t-1] + in * eff_i - out / eff_o')
    for (stf, sit, sto, com) in m_parameters['sto_tuples']:
        m.add_constraints(m.variables['e_sto_in'].loc[:, (stf, sit, sto, com)] -
                m_parameters['dt'] * m_parameters['cap_sto_p'][stf, sit, sto, com],
                "<=", 0, name="res_storage_input_by_power"+str((stf, sit, sto, com)))
        # m.res_storage_input_by_power = pyomo.Constraint(
        #     m.tm, m.sto_tuples,
        #     rule=res_storage_input_by_power_rule,
        #     doc='storage input <= storage power')
        m.add_constraints(m.variables['e_sto_out'].loc[:, (stf, sit, sto, com)] -
                m_parameters['dt'] * m_parameters['cap_sto_p'][stf, sit, sto, com],
                "<=", 0, name="res_storage_output_by_power"+str((stf, sit, sto, com)))
        # m.res_storage_output_by_power = pyomo.Constraint(
        #     m.tm, m.sto_tuples,
        #     rule=res_storage_output_by_power_rule,
        #     doc='storage output <= storage power')
        m.add_constraints(m.variables['e_sto_con'].loc[:, (stf, sit, sto, com)] -
                m_parameters['cap_sto_c'][stf, sit, sto, com],
                "<=", 0, name="res_storage_state_by_capacity"+str((stf, sit, sto, com)))
        # m.res_storage_state_by_capacity = pyomo.Constraint(
        #     m.t, m.sto_tuples,
        #     rule=res_storage_state_by_capacity_rule,
        #     doc='storage content <= storage capacity')
        m.add_constraints(m_parameters['cap_sto_p'][stf, sit, sto, com], "=>",
                m_parameters['storage_dict']['cap-lo-p'][stf, sit, sto, com],
                name="res_storage_power_low"+str((stf, sit, sto, com)))
        m.add_constraints(m_parameters['cap_sto_p'][stf, sit, sto, com], "<=",
                m_parameters['storage_dict']['cap-up-p'][stf, sit, sto, com],
                name="res_storage_power_high"+str((stf, sit, sto, com)))
        # m.res_storage_power = pyomo.Constraint(
        #     m.sto_tuples,
        #     rule=res_storage_power_rule,
        #     doc='storage.cap-lo-p <= storage power <= storage.cap-up-p')
        m.add_constraints(m_parameters['cap_sto_c'][stf, sit, sto, com], "=>",
                m_parameters['storage_dict']['cap-lo-c'][stf, sit, sto, com],
                name="res_storage_capacity_low"+str((stf, sit, sto, com)))
        m.add_constraints(m_parameters['cap_sto_c'][stf, sit, sto, com], "<=",
                m_parameters['storage_dict']['cap-up-c'][stf, sit, sto, com],
                name="res_storage_capacity_high"+str((stf, sit, sto, com)))
        # m.res_storage_capacity = pyomo.Constraint(
        #     m.sto_tuples,
        #     rule=res_storage_capacity_rule,
        #     doc='storage.cap-lo-c <= storage capacity <= storage.cap-up-c')
        m.add_constraints(m.variables['e_sto_con'].loc[:,(stf, sit, sto, com)][1,].reset_coords(drop=True) -
                m.variables['e_sto_con'].loc[:,(stf, sit, sto, com)][-1,].reset_coords(drop=True), "<=",
                0, name="res_storage_state_cyclicity"+str((stf, sit, sto, com)))
        # m.res_storage_state_cyclicity = pyomo.Constraint(
        #     m.sto_tuples,
        #     rule=res_storage_state_cyclicity_rule,
        #     doc='storage content initial <= final, both variable')
    for (stf, sit, sto, com) in m_parameters['sto_init_bound_tuples']:
        m.add_constraints(m.variables['e_sto_con'].loc[:,(stf, sit, sto, com)][1,] -
            m_parameters['cap_sto_c'][stf, sit, sto, com] *
            m_parameters['storage_dict']['init'][(stf, sit, sto, com)], "=", 0,
            name="def_initial_storage_state"+str((stf, sit, sto, com)))
        # m.def_initial_storage_state = pyomo.Constraint(
        #     m.sto_init_bound_tuples,
        #     rule=def_initial_storage_state_rule,
        #     doc='storage content initial == and final >= storage.init * capacity')
    for (stf, sit, sto, com) in m_parameters['sto_ep_ratio_tuples']:
        m.add_constraints(m_parameters['cap_sto_c'][stf, sit, sto, com] -
                m_parameters['cap_sto_p'][stf, sit, sto, com] *
                m_parameters['storage_dict']['ep-ratio'][(stf, sit, sto, com)],
                "=", 0, name="def_storage_energy_power_ratio"+str((stf, sit, sto, com)))
        # m.def_storage_energy_power_ratio = pyomo.Constraint(
        #     m.sto_ep_ratio_tuples,
        #     rule=def_storage_energy_power_ratio_rule,
        #     doc='storage capacity = storage power * storage E2P ratio')

    return m


# constraints

# storage content in timestep [t] == storage content[t-1] * (1-discharge)
# + newly stored energy * input efficiency
# - retrieved energy / output efficiency
def def_storage_state_rule(m, t, stf, sit, sto, com):
    return (m.e_sto_con[t, stf, sit, sto, com] ==
            m.e_sto_con[t - 1, stf, sit, sto, com] *
            (1 - m.storage_dict['discharge']
             [(stf, sit, sto, com)]) ** m.dt.value +
            m.e_sto_in[t, stf, sit, sto, com] *
            m.storage_dict['eff-in'][(stf, sit, sto, com)] -
            m.e_sto_out[t, stf, sit, sto, com] /
            m.storage_dict['eff-out'][(stf, sit, sto, com)])


# storage capacity (for m.cap_sto_c expression)
def def_storage_capacity_rule(m, m_parameters, stf, sit, sto, com):
    if m_parameters['mode']['int']:
        if (sit, sto, com, stf) in m_parameters['inst_sto_tuples']:
            if (min(m_parameters['stf']), sit, sto, com) in m_parameters['sto_const_cap_c_dict']:
                cap_sto_c = (m_parameters['storage_dict']['inst-cap-c'][
                    (min(m_parameters['stf']), sit, sto, com)] * 
                    m.variables['1'])
            else:
                cap_sto_c = (
                    sum(m.variables['cap_sto_c_new'].loc[(stf_built, sit, sto, com),].reset_coords(drop=True)
                        for stf_built in m.stf
                        if (sit, sto, com, stf_built, stf) in
                        m.operational_sto_tuples) +
                    m.storage_dict['inst-cap-c'][(min(m_parameters['stf']), sit, sto, com)])
        else:
            cap_sto_c = (
                sum(m.variables['cap_sto_c_new'].loc[(stf_built, sit, sto, com),].reset_coords(drop=True)
                    for stf_built in m.stf
                    if (sit, sto, com, stf_built, stf) in
                    m_parameters['operational_sto_tuples']))
    else:
        if (stf, sit, sto, com) in m_parameters['sto_const_cap_c_dict']:
            cap_sto_c = (
                    m_parameters['storage_dict']['inst-cap-c'][
                        (stf, sit, sto, com)] *
                    m.variables['1'])
        else:
            cap_sto_c = (m.variables['cap_sto_c_new'].loc[(stf, sit, sto, com),].reset_coords(drop=True) +
                         m_parameters['storage_dict']['inst-cap-c'][(stf, sit, sto, com)] *
                         m.variables['1'])

    return cap_sto_c


# storage power (for m.cap_sto_p expression)
def def_storage_power_rule(m, m_parameters, stf, sit, sto, com):
    if m_parameters['mode']['int']:
        if (sit, sto, com, stf) in m_parameters['inst_sto_tuples']:
            if (min(m_parameters['stf']), sit, sto, com) in m_parameters['sto_const_cap_p_dict']:
                cap_sto_p = (m_parameters['storage_dict']['inst-cap-p'][
                    (min(m_parameters['stf']), sit, sto, com)] *
                    m.variables['1'])
            else:
                cap_sto_p = (
                    sum(m.variables['cap_sto_p_new'].loc[(stf_built, sit, sto, com),].reset_coords(drop=True)
                        for stf_built in m_parameters['stf']
                        if (sit, sto, com, stf_built, stf) in
                        m_parameters['operational_sto_tuples']) +
                    m_parameters['storage_dict']['inst-cap-p'][
                        (min(m_parameters['stf']), sit, sto, com)] *
                    m.variables['1'])
        else:
            cap_sto_p = (
                sum(m.variables['cap_sto_p_new'].loc[(stf_built, sit, sto, com),].reset_coords(drop=True)
                    for stf_built in m_parameters['stf']
                    if (sit, sto, com, stf_built, stf)
                    in m_parameters['operational_sto_tuples']))
    else:
        if (stf, sit, sto, com) in m_parameters['sto_const_cap_p_dict']:
            cap_sto_p = (m_parameters['storage_dict']['inst-cap-p'][(stf, sit, sto, com)] *
                    m.variables['1'])
        else:
            cap_sto_p = (m.variables['cap_sto_p_new'].loc[(stf, sit, sto, com),].reset_coords(drop=True) +
                         m_parameters['storage_dict']['inst-cap-p'][(stf, sit, sto, com)] *
                         m.variables['1'])

    return cap_sto_p


# storage input <= storage power
def res_storage_input_by_power_rule(m, t, stf, sit, sto, com):
    return (m.e_sto_in[t, stf, sit, sto, com] <= m.dt *
            m.cap_sto_p[stf, sit, sto, com])


# storage output <= storage power
def res_storage_output_by_power_rule(m, t, stf, sit, sto, co):
    return (m.e_sto_out[t, stf, sit, sto, co] <= m.dt *
            m.cap_sto_p[stf, sit, sto, co])


# storage content <= storage capacity
def res_storage_state_by_capacity_rule(m, t, stf, sit, sto, com):
    return (m.e_sto_con[t, stf, sit, sto, com] <=
            m.cap_sto_c[stf, sit, sto, com])


# lower bound <= storage power <= upper bound
def res_storage_power_rule(m, stf, sit, sto, com):
    return (m.storage_dict['cap-lo-p'][(stf, sit, sto, com)],
            m.cap_sto_p[stf, sit, sto, com],
            m.storage_dict['cap-up-p'][(stf, sit, sto, com)])


# lower bound <= storage capacity <= upper bound
def res_storage_capacity_rule(m, stf, sit, sto, com):
    return (m.storage_dict['cap-lo-c'][(stf, sit, sto, com)],
            m.cap_sto_c[stf, sit, sto, com],
            m.storage_dict['cap-up-c'][(stf, sit, sto, com)])


# initialization of storage content in first timestep t[1]
# forced minimun  storage content in final timestep t[len(m.t)]
# content[t=1] == storage capacity * fraction <= content[t=final]
def def_initial_storage_state_rule(m, stf, sit, sto, com):
    return (m.e_sto_con[m.t[1], stf, sit, sto, com] ==
            m.cap_sto_c[stf, sit, sto, com] *
            m.storage_dict['init'][(stf, sit, sto, com)])


def res_storage_state_cyclicity_rule(m, stf, sit, sto, com):
    return (m.e_sto_con[m.t[1], stf, sit, sto, com] <=
            m.e_sto_con[m.t[len(m.t)], stf, sit, sto, com])


def def_storage_energy_power_ratio_rule(m, stf, sit, sto, com):
    return (m.cap_sto_c[stf, sit, sto, com] == m.cap_sto_p[stf, sit, sto, com] *
            m.storage_dict['ep-ratio'][(stf, sit, sto, com)])


# storage balance
def storage_balance(m, m_parameters, stf, sit, com):
    """called in commodity balance
    For a given commodity co , calculate the balance of
    storage input and output at all timesteps"""

   # usage as input for storage increases consumption
   # output from storage decreases consumption
    return (m.variables['e_sto_in'].loc[:,[(stframe, site, storage, com)
        for stframe, site, storage, commodity in m_parameters['sto_tuples']
        if site == sit and stframe == stf and commodity == com]].sum('sto_tuples') -
        m.variables['e_sto_out'].loc[:,[(stframe, site, storage, com)
        for stframe, site, storage, commodity in m_parameters['sto_tuples']
        if site == sit and stframe == stf and commodity == com]].sum('sto_tuples')) 


# storage costs
def storage_cost(m, m_parameters, cost_type):
    """returns storage cost function for the different cost types"""
    if cost_type == 'Invest':
        cost = sum(m.variables['cap_sto_p_new'].loc[s,].reset_coords(drop=True) *
                   m_parameters['storage_dict']['inv-cost-p'][s] *
                   m_parameters['storage_dict']['invcost-factor'][s] +
                   m.variables['cap_sto_c_new'].loc[s,].reset_coords(drop=True) *
                   m_parameters['storage_dict']['inv-cost-c'][s] *
                   m_parameters['storage_dict']['invcost-factor'][s]
                   for s in m_parameters['sto_tuples'])
        if m_parameters['mode']['int']:
            cost -= sum(m.variables['cap_sto_p_new'].loc[s,].reset_coords(drop=True) *
                        m_parameters['storage_dict']['inv-cost-p'][s] *
                        m_parameters['storage_dict']['overpay-factor'][s] +
                        m.variables['cap_sto_c_new'].loc[s,].reset_coords(drop=True) *
                        m_parameters['storage_dict']['inv-cost-c'][s] *
                        m_parameters['storage_dict']['overpay-factor'][s]
                        for s in m_parameters['sto_tuples'])
        return cost
    elif cost_type == 'Fixed':
        return sum((m_parameters['cap_sto_p'][s] *
                    m_parameters['storage_dict']['fix-cost-p'][s] +
                    m_parameters['cap_sto_c'][s] *
                    m_parameters['storage_dict']['fix-cost-c'][s]) *
                   m_parameters['storage_dict']['cost_factor'][s]
                   for s in m_parameters['sto_tuples'])
    elif cost_type == 'Variable':
        return sum(m.variables['e_sto_con'].loc[:, s].sum() *
                   m_parameters['weight'] *
                   m_parameters['storage_dict']['var-cost-c'][s] *
                   m_parameters['storage_dict']['cost_factor'][s] +
                   (m.variables['e_sto_in'].loc[:, s].sum() +
                       m.variables['e_sto_out'].loc[:, s].sum()) *
                   m_parameters['weight'] *
                   m_parameters['storage_dict']['var-cost-p'][s] *
                   m_parameters['storage_dict']['cost_factor'][s]
                   for s in m_parameters['sto_tuples'])


def op_sto_tuples(sto_tuple, m_parameters):
    """ s.a. op_pro_tuples
    """
    op_sto = []
    sorted_stf = sorted(list(m_parameters['stf']))

    for (stf, sit, sto, com) in sto_tuple:
        for stf_later in sorted_stf:
            index_helper = sorted_stf.index(stf_later)
            if stf_later == max(sorted_stf):
                if (stf_later +
                    m_parameters['global_prop_dict']['value'][(max(sorted_stf), 'Weight')] -
                    1 <= stf +
                        m_parameters['storage_dict']['depreciation'][(stf, sit, sto, com)]):
                    op_sto.append((sit, sto, com, stf, stf_later))
            elif (sorted_stf[index_helper + 1] <=
                  stf +
                  m_parameters['storage_dict']['depreciation'][(stf, sit, sto, com)] and
                  stf <= stf_later):
                op_sto.append((sit, sto, com, stf, stf_later))
            else:
                pass

    return op_sto


def inst_sto_tuples(m_parameters):
    """ s.a. inst_pro_tuples
    """
    inst_sto = []
    sorted_stf = sorted(list(m_parameters['stf']))

    for (stf, sit, sto, com) in m_parameters['inst_sto'].index:
        for stf_later in sorted_stf:
            index_helper = sorted_stf.index(stf_later)
            if stf_later == max(m_parameters['stf']):
                if (stf_later +
                    m_parameters['global_prop_dict']['value'][(max(sorted_stf), 'Weight')] -
                    1 < min(m_parameters['stf']) +
                        m_parameters['storage_dict']['lifetime'][(stf, sit, sto, com)]):
                    inst_sto.append((sit, sto, com, stf_later))
            elif (sorted_stf[index_helper + 1] <=
                  min(m_parameters['stf']) +
                  m_parameters['storage_dict']['lifetime'][
                                (stf, sit, sto, com)]):
                inst_sto.append((sit, sto, com, stf_later))

    return inst_sto
