import math
import pyomo.core as pyomo
import linopy
import pandas as pd
from .modelhelper import commodity_subset


def add_buy_sell_price(m, m_parameters):

    # Sets
    m_parameters['com_sell'] = pd.Index(commodity_subset(m_parameters['com_tuples'], 'Sell'),
            name="com_sell")
    # m.com_sell = pyomo.Set(
    #     within=m.com,
    #     initialize=commodity_subset(m.com_tuples, 'Sell'),
    #     doc='Commodities that can be sold')
    m_parameters['com_buy'] = pd.Index(commodity_subset(m_parameters['com_tuples'], 'Buy'),
            name="com_buy")
    # m.com_buy = pyomo.Set(
    #     within=m.com,
    #     initialize=commodity_subset(m.com_tuples, 'Buy'),
    #     doc='Commodities that can be purchased')

    # Variables
    m.add_variables(name="e_co_sell", coords=[m_parameters['tm'],
        m_parameters['com_tuples']], lower=0)
    # m.e_co_sell = pyomo.Var(
    #     m.tm, m.com_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='Use of sell commodity source (MW) per timestep')
    m.add_variables(name="e_co_buy", coords=[m_parameters['tm'],
        m_parameters['com_tuples']], lower=0)
    # m.e_co_buy = pyomo.Var(
    #     m.tm, m.com_tuples,
    #     within=pyomo.NonNegativeReals,
    #     doc='Use of buy commodity source (MW) per timestep')

    # Rules
    for (stf, sit, com, com_type) in m_parameters['com_tuples']:
        if com in m_parameters['com_sell']:
            m.add_constraints(m.variables['e_co_sell'].loc[:,(stf, sit, com, com_type)]
                    .reset_coords(drop=True),
                    "<=", m_parameters['dt'] *
                    m_parameters['commodity_dict']['maxperhour']
                    [(stf, sit, com, com_type)],
                    name="res_sell_step"+str((stf, sit, com, com_type)))
            # m.res_sell_step = pyomo.Constraint(
            #     m.tm, m.com_tuples,
            #     rule=res_sell_step_rule,
            #     doc='sell commodity output per step <= commodity.maxperstep')

            # calculate total sale of commodity com
            total_consumption = m.variables['e_co_sell'].loc[:,(stf,sit,com,com_type)].sum()
            total_consumption *= m_parameters['weight']
            m.add_constraints(total_consumption, "<=",
                    m_parameters['commodity_dict']['max']
                    [(stf, sit, com, com_type)],
                    name="res_sell_total"+str((stf, sit, com, com_type)))
            # m.res_sell_total = pyomo.Constraint(
            #     m.com_tuples,
            #     rule=res_sell_total_rule,
            #     doc='total sell commodity output <= commodity.max')

        if com in m_parameters['com_buy']:
            m.add_constraints(m.variables['e_co_buy'].loc[:,(stf, sit, com, com_type)]
                    .reset_coords(drop=True),
                    "<=", m_parameters['dt'] *
                    m_parameters['commodity_dict']['maxperhour']
                    [(stf, sit, com, com_type)],
                    name="res_buy_step"+str((stf, sit, com, com_type)))
            # m.res_buy_step = pyomo.Constraint(
            #     m.tm, m.com_tuples,
            #     rule=res_buy_step_rule,
            #     doc='buy commodity output per step <= commodity.maxperstep')
            # calculate total sale of commodity com
            total_consumption = m.variables['e_co_buy'].loc[:,(stf,sit,com,com_type)].sum()
            total_consumption *= m_parameters['weight']
            m.add_constraints(total_consumption, "<=",
                    m_parameters['commodity_dict']['max']
                    [(stf, sit, com, com_type)],
                    name="res_sell_total"+str((stf, sit, com, com_type)))
            # m.res_buy_total = pyomo.Constraint(
            #     m.com_tuples,
            #     rule=res_buy_total_rule,
            #     doc='total buy commodity output <= commodity.max')

    for (stf, sit_in, pro_in, coin) in m_parameters['pro_input_tuples']:
        if coin in m_parameters['com_buy']:
            sell_pro = search_sell_buy_tuple(m_parameters, stf, sit_in, pro_in, coin)
            if sell_pro is None:
                continue
            else:
                m.add_constraints(m_parameters['cap_pro'][stf, sit_in, pro_in] -
                        m_parameters['cap_pro'][stf, sit_in, sell_pro], "=", 0,
                        name="res_sell_buy_symmetry"+str((stf, sit_in, sell_pro)))
        # m.res_sell_buy_symmetry = pyomo.Constraint(
        #     m.pro_input_tuples,
        #     rule=res_sell_buy_symmetry_rule,
        #     doc='power connection capacity must be symmetric in both directions')

    return m


# constraints

# limit sell commodity use per time step
def res_sell_step_rule(m, tm, stf, sit, com, com_type):
    if com not in m.com_sell:
        return pyomo.Constraint.Skip
    else:
        return (m.e_co_sell[tm, stf, sit, com, com_type] <=
                m.dt * m.commodity_dict['maxperhour']
                [(stf, sit, com, com_type)])


# limit sell commodity use in total (scaled to annual consumption, thanks
# to m.weight)
def res_sell_total_rule(m, stf, sit, com, com_type):
    if com not in m.com_sell:
        return pyomo.Constraint.Skip
    else:
        # calculate total sale of commodity com
        total_consumption = 0
        for tm in m.tm:
            total_consumption += (
                m.e_co_sell[tm, stf, sit, com, com_type])
        total_consumption *= m.weight
        return (total_consumption <=
                m.commodity_dict['max'][(stf, sit, com, com_type)])


# limit buy commodity use per time step
def res_buy_step_rule(m, tm, stf, sit, com, com_type):
    if com not in m.com_buy:
        return pyomo.Constraint.Skip
    else:
        return (m.e_co_buy[tm, stf, sit, com, com_type] <=
                m.dt * m.commodity_dict['maxperhour']
                [(stf, sit, com, com_type)])


# limit buy commodity use in total (scaled to annual consumption, thanks
# to m.weight)
def res_buy_total_rule(m, stf, sit, com, com_type):
    if com not in m.com_buy:
        return pyomo.Constraint.Skip
    else:
        # calculate total sale of commodity com
        total_consumption = 0
        for tm in m.tm:
            total_consumption += (
                m.e_co_buy[tm, stf, sit, com, com_type])
        total_consumption *= m.weight
        return (total_consumption <=
                m.commodity_dict['max'][(stf, sit, com, com_type)])


# power connection capacity: Sell == Buy
def res_sell_buy_symmetry_rule(m, stf, sit_in, pro_in, coin):
    # constraint only for sell and buy processes
    # and the processes must be in the same site
    if coin in m.com_buy:
        sell_pro = search_sell_buy_tuple(m, stf, sit_in, pro_in, coin)
        if sell_pro is None:
            return pyomo.Constraint.Skip
        else:
            return (m.cap_pro[stf, sit_in, pro_in] ==
                    m.cap_pro[stf, sit_in, sell_pro])
    else:
        return pyomo.Constraint.Skip


def search_sell_buy_tuple(m_parameters, stf, sit_in, pro_in, coin):
    """ Return the equivalent sell-process for a given buy-process.
    Args:
        m: the parameters of the model
        sit_in: a site
        pro_in: a process
        co_in: a commodity
    Returns:
        a process
    """
    pro_output_tuples = [x for x in m_parameters['pro_output_tuples'] if x[1] == sit_in]
    pro_input_tuples = [x for x in m_parameters['pro_input_tuples'] if x[1] == sit_in]
    # search the output commodities for the "buy" process
    # buy_out = (stf, site, output_commodity)
    buy_out = set([(x[0], x[1], x[3])
                   for x in pro_output_tuples
                   if x[2] == pro_in])
    # search the sell process for the output_commodity from the buy process
    sell_output_tuple = ([x
                          for x in pro_output_tuples
                          if x[3] in m_parameters['com_sell']])
    for k in range(len(sell_output_tuple)):
        sell_pro = sell_output_tuple[k][2]
        sell_in = set([(x[0], x[1], x[3])
                       for x in pro_input_tuples
                       if x[2] == sell_pro])
        # check: buy - commodity == commodity - sell; for a site
        if not(sell_in.isdisjoint(buy_out)):
            return sell_pro
    return None


def bsp_surplus(m, m_parameters, stf, sit, com, com_type):

    power_surplus = linopy.LinearExpression()

    # if com is a sell commodity, the commodity source term e_co_sell
    # can supply a possibly positive power_surplus
    if com in m_parameters['com_sell']:
        power_surplus -= m.variables['e_co_sell'].loc[:,
                (stf, sit, com, com_type)].reset_coords(drop=True)

    # if com is a buy commodity, the commodity source term e_co_buy
    # can supply a possibly negative power_surplus
    if com in m_parameters['com_buy']:
        power_surplus += m.variables['e_co_buy'].loc[:,
                (stf, sit, com, com_type)].reset_coords(drop=True)

    return power_surplus


def revenue_costs(m, m_parameters):
    sell_tuples = commodity_subset(m_parameters['com_tuples'],
            m_parameters['com_sell'])
    try:
        return -sum(
            m.variables['e_co_sell'].loc[tm, c].reset_coords(drop=True) *
            m_parameters['buy_sell_price_dict'][c[2]][(c[0], tm)] *
            m_parameters['weight'] *
            m_parameters['commodity_dict']['price'][c] *
            m_parameters['commodity_dict']['cost_factor'][c]
            for tm in m_parameters['tm']
            for c in sell_tuples)
    except KeyError:
        return -sum(
            m.variables['e_co_sell'].loc[tm, c].reset_coords(drop=True) *
            m_parameters['buy_sell_price_dict'][c[2], ][(c[0], tm)] * 
            m_parameters['weight'] *
            m_parameters['commodity_dict']['price'][c] *
            m_parameters['commodity_dict']['cost_factor'][c]
            for tm in m_parameters['tm']
            for c in sell_tuples)


def purchase_costs(m, m_parameters):
    buy_tuples = commodity_subset(m_parameters['com_tuples'], m_parameters['com_buy'])
    try:
        return sum(
            m.variables['e_co_buy'].loc[tm, c].reset_coords(drop=True) *
            m_parameters['buy_sell_price_dict'][c[2]][(c[0], tm)] * 
            m_parameters['weight'] *
            m_parameters['commodity_dict']['price'][c] *
            m_parameters['commodity_dict']['cost_factor'][c]
            for tm in m_parameters['tm']
            for c in buy_tuples)
    except KeyError:
        return sum(
            m.variables['e_co_buy'].loc[tm, c].reset_coords(drop=True) *
            m_parameters['buy_sell_price_dict'][c[2], ][(c[0], tm)] * 
            m_parameters['weight'] *
            m_parameters['commodity_dict']['price'][c] *
            m_parameters['commodity_dict']['cost_factor'][c]
            for tm in m_parameters['tm']
            for c in buy_tuples)
