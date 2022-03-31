import math
import pyomo.core as pyomo
import pandas as pd
import linopy

def e_tra_domain_rule(m_parameters, tm, stf, sin, sout, tra, com):
    # assigning e_tra_in and e_tra_out variable domains for transport and DCPF
    if (stf, sin, sout, tra, com) in m_parameters['tra_tuples_dc']:
        return pyomo.Reals
    elif (stf, sin, sout, tra, com) in m_parameters['tra_tuples_tp']:
        return pyomo.NonNegativeReals


def remove_duplicate_transmission(transmission_keys):
    # removing duplicate transmissions for DCPF
    tra_tuple_list = list(transmission_keys)
    i = 0
    while i < len(tra_tuple_list):
        for k in range(len(tra_tuple_list)):
            if (tra_tuple_list[i][1] == tra_tuple_list[k][2] and
                    tra_tuple_list[i][2] == tra_tuple_list[k][1] and
                    tra_tuple_list[i][0] == tra_tuple_list[k][0] and
                    tra_tuple_list[i][3] == tra_tuple_list[k][3]):
                del tra_tuple_list[i]
                i -= 1
                break
        i += 1
    return set(tra_tuple_list)


def add_transmission(m, m_parameters):

    # tranmission (e.g. hvac, hvdc, pipeline...)
    indexlist = set()
    for key in m_parameters['transmission_dict']["eff"]:
        indexlist.add(tuple(key)[3])
    m_parameters['tra'] = pd.Index(indexlist, name="tra")
    # doc='Set of transmission technologies')

    # transmission tuples
    m_parameters['tra_tuples'] = pd.Index(m_parameters['transmission_dict']["eff"].keys(),
            name="tra_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.sit * m.tra * m.com,
    # doc='Combinations of possible transmissions, e.g. '
    #     '(2020,South,Mid,hvac,Elec)')

    if m_parameters['mode']['int']:
        m_parameters['operational_tra_tuples'] = pd.Index([(sit, sit_, tra, com, stf, stf_later)
                                    for (sit, sit_, tra, com, stf, stf_later)
                                    in op_tra_tuples(m_parameters['tra_tuples'], m_parameters)],
                                    name="operational_tra_tuples", tupleize_cols=False)
        # within=m.sit * m.sit * m.tra * m.com * m.stf * m.stf,
        # doc='Transmissions that are still operational through stf_later'
        #     '(and the relevant years following), if built in stf'
        #     'in stf.')
        m_parameters['inst_tra_tuples'] = pd.Index([(sit, sit_, tra, com, stf)
                             for (sit, sit_, tra, com, stf)
                             in inst_tra_tuples(m_parameters)],
                             name="inst_tra_tuples", tupleize_cols=False)
        # within=m.sit * m.sit * m.tra * m.com * m.stf,
        # doc='Installed transmissions that are still operational'
        #     'through stf')

    # Variables
    m.add_variables(lower=0, coords=[m_parameters['tra_tuples']], name="cap_tra_new")
    # doc='New transmission capacity (MW)')

    # transmission capacity as expression object
    # m.cap_tra = pyomo.Expression(
    #     m.tra_tuples,
    #     rule=def_transmission_capacity_rule,
    #     doc='total transmission capacity')
    # TODO : Try to have this as a single linopy.LinearExpression, instead of
    # a dict of linopy.LinearExpression
    m_parameters['cap_tra'] = {(stf, sin, sout, tra, com): 
            def_transmission_capacity_rule(m, m_parameters, stf, sin, sout, tra, com)
            for (stf, sin, sout, tra, com) in m_parameters['tra_tuples']}

    m.add_variables(lower=0, coords=[m_parameters['tm'], m_parameters['tra_tuples']], name="e_tra_in")
    # doc='Power flow into transmission line (MW) per timestep')
    m.add_variables(lower=0, coords=[m_parameters['tm'], m_parameters['tra_tuples']], name="e_tra_out")
    # doc='Power flow out of transmission line (MW) per timestep')

    # transmission
    m.add_constraints(
            m.variables['e_tra_out'] - 
            m.variables['e_tra_in'] *
            pd.Series(m_parameters['transmission_dict']['eff'], name="tra_tuples") == 0,
            name="def_transmission_output")
    # m.def_transmission_output = pyomo.Constraint(
    #     m.tm, m.tra_tuples,
    #     rule=def_transmission_output_rule,
    #     doc='transmission output = transmission input * efficiency')
    for (stf, sin, sout, tra, com) in m_parameters['tra_tuples']:
        m.add_constraints(
                m.variables['e_tra_in'].loc[:,(stf, sin, sout, tra, com)] -
                m_parameters['dt'] *
                m_parameters['cap_tra'][(stf, sin, sout, tra, com)] <= 0,
                name="res_transmission_input_by_capacity"+str((stf, sin, sout, tra, com)))
        # m.res_transmission_input_by_capacity = pyomo.Constraint(
        #     m.tm, m.tra_tuples,
        #     rule=res_transmission_input_by_capacity_rule,
        #     doc='transmission input <= total transmission capacity')
        m.add_constraints(
                m_parameters['cap_tra'][(stf, sin, sout, tra, com)], ">=",
                m_parameters['transmission_dict']['cap-lo'][
                    (stf, sin, sout, tra, com)],
                name="res_transmission_capacity_low"+str(
                    (stf, sin, sout, tra, com)))
        m.add_constraints(
                m_parameters['cap_tra'][(stf, sin, sout, tra, com)], "<=", 
                m_parameters['transmission_dict']['cap-up'][(stf, sin, sout, tra, com)],
                name="res_transmission_capacity_high"+str((stf, sin, sout, tra, com)))
        # m.res_transmission_capacity = pyomo.Constraint(
        #     m.tra_tuples,
        #     rule=res_transmission_capacity_rule,
        #     doc='transmission.cap-lo <= total transmission capacity <= '
        #         'transmission.cap-up')
    # TODO: Find a way to exchange sin and sout coordinates using a linear
    # expression (Some sort of transposition, is it considered a linear
    # expression by linopy?). There is a method transpose() that returns a
    # variable with transposed coordinates. That would require the tuples to be
    # considered as multiple dimensions.a That pose other problems, though.
    for (stf, sin, sout, tra, com) in m_parameters['tra_tuples']:
        m.add_constraints(
            m_parameters['cap_tra'][stf, sin, sout, tra, com] -
            m_parameters['cap_tra'][stf, sout, sin, tra, com],
            "=",
            0,
            name="res_transmission_symmetry"+str((stf, sin, sout, tra, com)))
    # m.res_transmission_symmetry = pyomo.Constraint(
    #     m.tra_tuples,
    #     rule=res_transmission_symmetry_rule,
    #     doc='total transmission capacity must be symmetric in both directions')

    return m

# adds the transmission features to model with DCPF model features
def add_transmission_dc(m, m_parameters):
    # defining transmission tuple sets for transport and DCPF model separately
    tra_tuples = set()
    tra_tuples_dc = set()
    for key in m_parameters['transmission_dict']['reactance']:
        tra_tuples.add(tuple(key))
    for key in m_parameters['transmission_dc_dict']['reactance']:
        tra_tuples_dc.add(tuple(key))
    tra_tuples_tp = tra_tuples - tra_tuples_dc
    tra_tuples_dc = remove_duplicate_transmission(tra_tuples_dc)
    tra_tuples = tra_tuples_dc | tra_tuples_tp

    # tranmission (e.g. hvac, hvdc, pipeline...)
    indexlist = set()
    for key in m_parameters['transmission_dict']["eff"]:
        indexlist.add(tuple(key)[3])
    m_parameters['tra'] = pd.Index(indexlist, name="tra")
    # m.tra = pyomo.Set(
    #     initialize=indexlist,
    #     doc='Set of transmission technologies')

    # Transport and DCPF transmission tuples
    m_parameters['tra_tuples'] = pd.Index(tra_tuples, name="tra_tuples", 
            tupleize_cols=False)
    # m.tra_tuples = pyomo.Set(
    #     within=m.stf * m.sit * m.sit * m.tra * m.com,
    #     initialize=tuple(tra_tuples),
    #     doc='Combinations of possible transmissions,'
    #         'without duplicate dc transmissions'
    #         ' e.g. (2020,South,Mid,hvac,Elec)')

    # DCPF transmission tuples
    m_parameters['tra_tuples_dc'] = pd.Index(tra_tuples_dc, name="tra_tuples_dc",
            tupleize_cols=False)
    # m.tra_tuples_dc = pyomo.Set(
    #     within=m.stf * m.sit * m.sit * m.tra * m.com,
    #     initialize=tuple(tra_tuples_dc),
    #     doc='Combinations of possible bidirectional dc'
    #         'transmissions, e.g. (2020,South,Mid,hvac,Elec)')

    # Transport transmission tuples
    m_parameters['tra_tuples_tp'] = pd.Index(tra_tuples_tp, name="tra_tuples_tp",
            tupleize_cols=False)
    # m.tra_tuples_tp = pyomo.Set(
    #     within=m.stf * m.sit * m.sit * m.tra * m.com,
    #     initialize=tuple(tra_tuples_tp),
    #     doc='Combinations of possible transport transmissions,'
    #         'e.g. (2020,South,Mid,hvac,Elec)')

    if m_parameters['mode']['int']:
        m_parameters['operational_tra_tuples'] = pd.Index([(sit, sit_, tra, com, stf, stf_later)
                                    for (sit, sit_, tra, com, stf, stf_later)
                                    in op_tra_tuples(m_parameters['tra_tuples'], m_parameters)],
                                    name="operational_tra_tuples", tupleize_cols=False)
        # m.operational_tra_tuples = pyomo.Set(
        #     within=m.sit * m.sit * m.tra * m.com * m.stf * m.stf,
        #     initialize=[(sit, sit_, tra, com, stf, stf_later)
        #                 for (sit, sit_, tra, com, stf, stf_later)
        #                 in op_tra_tuples(m.tra_tuples, m)],
        #     doc='Transmissions that are still operational through stf_later'
        #         '(and the relevant years following), if built in stf'
        #         'in stf.')
        m_parameters['inst_tra_tuples'] = pd.Index([(sit, sit_, tra, com, stf)
                             for (sit, sit_, tra, com, stf)
                             in inst_tra_tuples(m)],
                             name="inst_tra_tuples", tupleize_cols=False)
        # m.inst_tra_tuples = pyomo.Set(
        #     within=m.sit * m.sit * m.tra * m.com * m.stf,
        #     initialize=[(sit, sit_, tra, com, stf)
        #                 for (sit, sit_, tra, com, stf)
        #                 in inst_tra_tuples(m)],
        #     doc='Installed transmissions that are still operational'
        #         'through stf')

    # Variables
    m.add_variables(lower=0, coords=[m_parameters['tra_tuples']],
            name="cap_tra_new")
    # m.cap_tra_new = pyomo.Var(
    #     m.tra_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='New transmission capacity (MW)')

    # transmission capacity as expression object
    m_parameters['cap_tra'] = {(stf, sin, sout, tra, com): 
            def_transmission_capacity_rule(m, m_parameters, stf, sin, sout, tra, com)
            for (stf, sin, sout, tra, com) in m_parameters['tra_tuples']}
    # m.cap_tra = pyomo.Expression(
    #     m.tra_tuples,
    #     rule=def_transmission_capacity_rule,
    #     doc='total transmission capacity')

    m.add_variables(lower=0, coords=[m_parameters['tm'], m_parameters['tra_tuples_dc']],
            name="e_tra_abs")
    # m.e_tra_abs = pyomo.Var(
    #     m.tm, m.tra_tuples_dc,
    #     within=pyomo.NonNegativeReals,
    #     doc='Absolute power flow on transmission line (MW) per timestep')
    #TODO : variable domain depending on whether the coordinate is in tra_tupples_dc or tra_tuples_tp
    # Maybe with constraints ?
    m.add_variables(coords=[m_parameters['tm'], m_parameters['tra_tuples']], name="e_tra_in")
    # m.e_tra_in = pyomo.Var(
    #     m.tm, m.tra_tuples,
    #     within=e_tra_domain_rule,
    #     doc='Power flow into transmission line (MW) per timestep')
    m.add_variables(coords=[m_parameters['tm'], m_parameters['tra_tuples']], name="e_tra_out")
    # m.e_tra_out = pyomo.Var(
    #     m.tm, m.tra_tuples,
    #     within=e_tra_domain_rule,
    #     doc='Power flow out of transmission line (MW) per timestep')

    m.add_variables(coords=[m_parameters['tm'], m_parameters['stf'],
        m_parameters['sit']], name="voltage_angle")
    # m.voltage_angle = pyomo.Var(
    #     m.tm, m.stf, m.sit,
    #     within=pyomo.Reals,
    #     doc='Voltage angle of a site')

    # transmission
    m.add_constraints(
            m.variables['e_tra_out'] -
            m.variables['e_tra_in'] *
            pd.Series(m_parameters['transmission_dict']['eff'], name="tra_tuples") == 0,
            name="def_transmission_output")
    # m.def_transmission_output = pyomo.Constraint(
    #     m.tm, m.tra_tuples,
    #     rule=def_transmission_output_rule,
    #     doc='transmission output = transmission input * efficiency')
    for (stf, sin, sout, tra, com) in m_parameters['tra_tuples_dc']:
        m.add_constraints(m.variables['e_tra_in'].loc[:, (stf, sin, sout, tra, com)] + 
            (m.variables['voltage_angle'].loc[:, stf, sin].reset_coords(drop=True) -
                m.variables['voltage_angle'].loc[:, stf, sout]).reset_coords(drop=True) / 57.2958 * 
                (-1 / m_parameters['transmission_dict']['reactance'][(stf, sin, sout, tra, com)])
                * m_parameters['transmission_dict']['base_voltage'][(stf, sin, sout, tra, com)]
                * m_parameters['transmission_dict']['base_voltage'][(stf, sin, sout, tra, com)],
                "=", 0,
                name="def_dc_power_flow"+str((stf, sin, sout, tra, com)))
        # m.def_dc_power_flow = pyomo.Constraint(
        #     m.tm, m.tra_tuples_dc,
        #     rule=def_dc_power_flow_rule,
        #     doc='transmission output = (angle(in)-angle(out))/ 57.2958 '
        #         '* -1 *(-1/reactance) * (base voltage)^2')
        m.add_constraints(- m_parameters['transmission_dict']['difflimit'][
            (stf, sin, sout, tra, com)] <=
            (m.variables['voltage_angle'].loc[:, stf, sin].reset_coords(drop=True) -
                m.variables['voltage_angle'].loc[:, stf, sout].reset_coords(drop=True)),
            name="def_angle_limit_low"+str((stf, sin, sout, tra, com)))
        m.add_constraints((m.variables['voltage_angle'].loc[:, stf, sin].reset_coords(drop=True) -
            m.variables['voltage_angle'].loc[:, stf, sout].reset_coords(drop=True)) <=
            m_parameters['transmission_dict']['difflimit'][
                (stf, sin, sout, tra, com)],
                name="def_angle_limit_high"+str((stf, sin, sout, tra, com)))
        # m.def_angle_limit = pyomo.Constraint(
        #     m.tm, m.tra_tuples_dc,
        #     rule=def_angle_limit_rule,
        #     doc='-angle limit < angle(in) - angle(out) < angle limit')
        m.add_constraints(- m.variables['e_tra_in'].loc[:, (stf, sin, sout, tra, com)] -
                m_parameters['dt'] * m_parameters['cap_tra'][stf, sin, sout, tra, com], "<=", 0,
                name="res_transmission_dc_input_by_capacity"+str((stf, sin, sout, tra, com)))
        # m.res_transmission_dc_input_by_capacity = pyomo.Constraint(
        #     m.tm, m.tra_tuples_dc,
        #     rule=res_transmission_dc_input_by_capacity_rule,
        #     doc='-dcpf transmission input <= total transmission capacity')

    m.add_constraints(m.variables['e_tra_in'] -
            m.variables['e_tra_abs'] <= 0,
            name="e_tra_abs1")
    m.add_constraints(-m.variables['e_tra_in'] -
            m.variables['e_tra_abs'] <= 0,
            name="e_tra_abs2")
    # m.e_tra_abs1 = pyomo.Constraint(
    #     m.tm, m.tra_tuples_dc,
    #     rule=e_tra_abs_rule1,
    #     doc='transmission dc input <= absolute transmission dc input')
    # m.e_tra_abs2 = pyomo.Constraint(
    #     m.tm, m.tra_tuples_dc,
    #     rule=e_tra_abs_rule2,
    #     doc='-transmission dc input <= absolute transmission dc input')

    for (stf, sin, sout, tra, com) in m_parameters['tra_tuples']:
        m.add_constraints(m.variables['e_tra_in'].loc[:, (stf, sin, sout, tra, com)] -
            m_parameters['dt'] * m_parameters['cap_tra'][stf, sin, sout, tra, com],
            "<=", 0, name="res_transmission_input_by_capacity"+
            str((stf, sin, sout, tra, com)))
        # m.res_transmission_input_by_capacity = pyomo.Constraint(
        #     m.tm, m.tra_tuples,
        #     rule=res_transmission_input_by_capacity_rule,
        #     doc='transmission input <= total transmission capacity')
        m.add_constraints(
                m_parameters['cap_tra'][(stf, sin, sout, tra, com)], ">=",
                m_parameters['transmission_dict']['cap-lo'][
                    (stf, sin, sout, tra, com)],
                name="res_transmission_capacity_low"+str(
                    (stf, sin, sout, tra, com)))
        m.add_constraints(
                m_parameters['cap_tra'][(stf, sin, sout, tra, com)], "<=", 
                m_parameters['transmission_dict']['cap-up'][(stf, sin, sout, tra, com)],
                name="res_transmission_capacity_high"+str((stf, sin, sout, tra, com)))
        # m.res_transmission_capacity = pyomo.Constraint(
        #     m.tra_tuples,
        #     rule=res_transmission_capacity_rule,
        #     doc='transmission.cap-lo <= total transmission capacity <= '
        #         'transmission.cap-up')

    for (stf, sin, sout, tra, com) in m_parameters['tra_tuples_tp']:
        m.add_constraints(
            m_parameters['cap_tra'][stf, sin, sout, tra, com] -
            m_parameters['cap_tra'][stf, sout, sin, tra, com],
            "=", 0,
            name="res_transmission_symmetry"+str((stf, sin, sout, tra, com)))
        # m.res_transmission_symmetry = pyomo.Constraint(
        #     m.tra_tuples_tp,
        #     rule=res_transmission_symmetry_rule,
        #     doc='total transmission capacity must be symmetric in both directions')

    return m


# constraints

# transmission capacity (for m.cap_tra expression)
def def_transmission_capacity_rule(m, m_parameters, stf, sin, sout, tra, com):
    if m_parameters['mode']['int']:
        if (sin, sout, tra, com, stf) in m_parameters['inst_tra_tuples']:
            if (min(m_parameters['stf']), sin, sout, tra, com) in m_parameters['tra_const_cap_dict']:
                cap_tra = (m_parameters['transmission_dict']['inst-cap'][
                    (min(m_parameters['stf']), sin, sout, tra, com)] *
                    m.variables['1'])
            else:
                cap_tra = (
                        sum(m.variables['cap_tra_new'].loc[(stf_built, sin,
                            sout, tra, com),].reset_coords(drop=True)
                        for stf_built in m_parameters['stf']
                        if (sin, sout, tra, com, stf_built, stf) in
                        m_parameters['operational_tra_tuples']) +
                    m_parameters['transmission_dict']['inst-cap']
                    [(min(m.stf), sin, sout, tra, com)] *
                    m.variables['1'])
        else:
            cap_tra = (
                    sum(m.variables['cap_tra_new'].loc[(stf_built, sin, sout,
                        tra, com),].reset_coords(drop=True)
                    for stf_built in m_parameters['stf']
                    if (sin, sout, tra, com, stf_built, stf) in
                    m_parameters['operational_tra_tuples']))
    else:
        if (stf, sin, sout, tra, com) in m_parameters['tra_const_cap_dict']:
            cap_tra = (m_parameters['transmission_dict']
                                   ['inst-cap']
                                   [(stf, sin, sout, tra, com)] *
                       m.variables["1"])
        else:
            cap_tra = (m.variables['cap_tra_new'].loc[(stf, sin, sout, tra,
                       com),].reset_coords(drop=True) +
                       m_parameters['transmission_dict']['inst-cap'][
                           (stf, sin, sout, tra, com)] *
                       m.variables["1"])

    return cap_tra

# transmission output == transmission input * efficiency


def def_transmission_output_rule(m, tm, stf, sin, sout, tra, com):
    return (m.e_tra_out[tm, stf, sin, sout, tra, com] ==
            m.e_tra_in[tm, stf, sin, sout, tra, com] *
            m.transmission_dict['eff'][(stf, sin, sout, tra, com)])

# power flow rule for DCPF transmissions
def def_dc_power_flow_rule(m, tm, stf, sin, sout, tra, com):
    return (m.e_tra_in[tm, stf, sin, sout, tra, com] ==
            (m.voltage_angle[tm, stf, sin] - m.voltage_angle[tm, stf, sout]) / 57.2958 * -1 *
            (-1 / m.transmission_dict['reactance'][(stf, sin, sout, tra, com)])
            * m.transmission_dict['base_voltage'][(stf, sin, sout, tra, com)]
            * m.transmission_dict['base_voltage'][(stf, sin, sout, tra, com)])

# voltage angle difference rule for DCPF transmissions
def def_angle_limit_rule(m, tm, stf, sin, sout, tra, com):
    return (- m.transmission_dict['difflimit'][(stf, sin, sout, tra, com)],
            (m.voltage_angle[tm, stf, sin] - m.voltage_angle[tm, stf, sout]),
            m.transmission_dict['difflimit'][(stf, sin, sout, tra, com)])

# first rule for creating absolute transmission input
def e_tra_abs_rule1(m, tm, stf, sin, sout, tra, com):
    return (m.e_tra_in[tm, stf, sin, sout, tra, com] <=
            m.e_tra_abs[tm, stf, sin, sout, tra, com])

# second rule for creating absolute transmission input
def e_tra_abs_rule2(m, tm, stf, sin, sout, tra, com):
    return (-m.e_tra_in[tm, stf, sin, sout, tra, com] <=
            m.e_tra_abs[tm, stf, sin, sout, tra, com])


# transmission input <= transmission capacity
def res_transmission_input_by_capacity_rule(m, tm, stf, sin, sout, tra, com):
    return (m.e_tra_in[tm, stf, sin, sout, tra, com] <=
            m.dt * m.cap_tra[stf, sin, sout, tra, com])


# - dc transmission input <= transmission capacity
def res_transmission_dc_input_by_capacity_rule(m, tm, stf, sin, sout, tra, com):
    return (- m.e_tra_in[tm, stf, sin, sout, tra, com] <=
            m.dt * m.cap_tra[stf, sin, sout, tra, com])


# lower bound <= transmission capacity <= upper bound
def res_transmission_capacity_rule(m, stf, sin, sout, tra, com):
    return (m.transmission_dict['cap-lo'][(stf, sin, sout, tra, com)],
            m.cap_tra[stf, sin, sout, tra, com],
            m.transmission_dict['cap-up'][(stf, sin, sout, tra, com)])


# transmission capacity from A to B == transmission capacity from B to A
def res_transmission_symmetry_rule(m, stf, sin, sout, tra, com):
    return m.cap_tra[stf, sin, sout, tra, com] == (m.cap_tra
                                                   [stf, sout, sin, tra, com])


# transmission balance
def transmission_balance(m, m_parameters, stf, sit, com):
    """called in commodity balance
    For a given commodity co and timestep tm, calculate the balance of
    import and export """

    return (m.variables['e_tra_in'].loc[:,[(stframe, site_in, site_out,
                            transmission, com)
                # exports increase balance
                for stframe, site_in, site_out, transmission, commodity
                in m_parameters['tra_tuples']
                if (site_in == sit and stframe == stf and
                    commodity == com)]].sum('tra_tuples') -
                m.variables['e_tra_out'].loc[:,[(stframe, site_in, site_out,
                             transmission, com)
                # imports decrease balance
                for stframe, site_in, site_out, transmission, commodity
                in m_parameters['tra_tuples']
                if (site_out == sit and stframe == stf and
                    commodity == com)]].sum('tra_tuples'))


# transmission cost function
def transmission_cost(m, m_parameters, cost_type):
    """returns transmission cost function for the different cost types"""
    if cost_type == 'Invest':
        cost = linopy.expressions.merge([m.variables['cap_tra_new'].loc[t,].reset_coords(drop=True)
            .reset_coords(drop=True) *
            m_parameters['transmission_dict']['inv-cost'][t] *
            m_parameters['transmission_dict']['invcost-factor'][t]
            for t in m_parameters['tra_tuples']])
        if m_parameters['mode']['int']:
            cost = linopy.expressions.merge([cost] + [
                -m.variables['cap_tra_new'].loc[t,].reset_coords(drop=True) *
                m_parameters['transmission_dict']['inv-cost'][t] *
                m_parameters['transmission_dict']['overpay-factor'][t]
                for t in m_parameters['tra_tuples']])
        return cost
    elif cost_type == 'Fixed':
        return sum(m_parameters['cap_tra'][t] * m_parameters['transmission_dict']['fix-cost'][t] *
                   m_parameters['transmission_dict']['cost_factor'][t]
                   for t in m_parameters['tra_tuples'])
    elif cost_type == 'Variable':
        if m_parameters['mode']['dpf']:
            return sum(m.variables['e_tra_in'].loc[:, t].sum() *
                       m_parameters['weight'] *
                       m_parameters['transmission_dict']['var-cost'][t] *
                       m_parameters['transmission_dict']['cost_factor'][t]
                       for t in m_parameters['tra_tuples_tp']) + \
                   sum(m.variables['e_tra_abs'].loc[:, t].sum() *
                       m_parameters['weight'] *
                       m_parameters['transmission_dict']['var-cost'][t] *
                       m_parameters['transmission_dict']['cost_factor'][t]
                       for t in m_parameters['tra_tuples_dc'])
        else:
            return sum(m.variables['e_tra_in'].loc[:, t].sum() *
                       m_parameters['weight'] *
                       m_parameters['transmission_dict']['var-cost'][t] *
                       m_parameters['transmission_dict']['cost_factor'][t]
                       for t in m_parameters['tra_tuples'])


def op_tra_tuples(tra_tuple, m_parameters):
    """ s.a. op_pro_tuples
    """
    op_tra = []
    sorted_stf = sorted(list(m_parameters['stf']))

    for (stf, sit1, sit2, tra, com) in tra_tuple:
        for stf_later in sorted_stf:
            index_helper = sorted_stf.index(stf_later)
            if stf_later == max(sorted_stf):
                if (stf_later +
                    m_parameters['global_prop_dict']['value'][(max(sorted_stf), 'Weight')] -
                    1 <= stf + m_parameters['transmission_dict']['depreciation'][
                        (stf, sit1, sit2, tra, com)]):
                    op_tra.append((sit1, sit2, tra, com, stf, stf_later))
            elif (sorted_stf[index_helper + 1] <=
                  stf + m_parameters['transmission_dict']['depreciation'][
                      (stf, sit1, sit2, tra, com)] and stf <= stf_later):
                op_tra.append((sit1, sit2, tra, com, stf, stf_later))
            else:
                pass

    return op_tra


def inst_tra_tuples(m_parameters):
    """ s.a. inst_pro_tuples
    """
    inst_tra = []
    sorted_stf = sorted(list(m_parameters['stf']))

    for (stf, sit1, sit2, tra, com) in m_parameters['inst_tra'].index:
        for stf_later in sorted_stf:
            index_helper = sorted_stf.index(stf_later)
            if stf_later == max(m_parameters['stf']):
                if (stf_later +
                    m_parameters['global_prop_dict']['value'][(max(sorted_stf), 'Weight')] -
                    1 < min(m_parameters['stf']) + m_parameters['transmission_dict']['lifetime'][
                        (stf, sit1, sit2, tra, com)]):
                    inst_tra.append((sit1, sit2, tra, com, stf_later))
            elif (sorted_stf[index_helper + 1] <= min(m_parameters['stf']) +
                  m_parameters['transmission_dict']['lifetime'][
                      (stf, sit1, sit2, tra, com)]):
                inst_tra.append((sit1, sit2, tra, com, stf_later))

    return inst_tra
