import math
import pyomo.core as pyomo
import linopy
from datetime import datetime
from .features import *
from .input import *


def create_model(data, dt=1, timesteps=None, objective='cost',
                 dual=True):
    """Create a linopy Model urbs object from given input data.

    Args:
        - data: a dict of up to 12
        - dt: timestep duration in hours (default: 1)
        - timesteps: optional list of timesteps, default: demand timeseries
        - objective: Either "cost" or "CO2" for choice of objective function,
          default: "cost"
        - dual: set True to add dual variables to model output
          (marginally slower), default: True

    Returns:
        a linopy Model object
    """

    # Optional
    if not timesteps:
        timesteps = data['demand'].index.tolist()
    m, m_parameters = linopy_model_prep(data, timesteps)  # preparing linopy model
    m_parameters['name'] = 'urbs'
    m_parameters['created'] = datetime.now().strftime('%Y%m%dT%H%M')
    m_parameters['_data'] = data

    # Parameters
    # Using Pyomo, the parameters were stored in pyomo.Param and pyomo.Set.
    # There are no such class in linopy. The parameters are directly stored
    # as parameters of the model object.
    # Warning : There is no validation like the one provided by the pyomo
    # classes -> when needed, the validation should be done manually

    # weight = length of year (hours) / length of simulation (hours)
    # weight scales costs and emissions from length of simulation to a full
    # year, making comparisons among cost types (invest is annualized, fixed
    # costs are annual by default, variable costs are scaled by weight) and
    # among different simulation durations meaningful.
    m_parameters['weight'] = (float(8760) /
            ((len(m_parameters['timesteps']) - 1) * dt))
    # doc='Pre-factor for variable costs and emissions for an annual result'

    # dt = spacing between timesteps. Required for storage equation that
    # converts between energy (storage content, e_sto_con) and power (all other
    # quantities that start with "e_")
    m_parameters['dt'] = dt
    # doc='Time step duration (in hours), default: 1'

    # import objective function information
    m_parameters['obj'] = objective
    # doc='Specification of minimized quantity, default: "cost"'

    # Sets
    # ====
    # Syntax: m.{name} = Set({domain}, initialize={values})
    # where name: set name
    #       domain: set domain for tuple sets, a cartesian set product
    #       values: set values, a list or array of element tuples

    # With linopy, there is no pyomo.Set class. The set parameters are stored
    # in native python iterable.
    # TODO: Check that the types are suitable (some lists should maybe be sets ?)

    # generate ordered time step sets
    m_parameters['t'] = m_parameters['timesteps']
    # doc='Set of timesteps'

    # modelled (i.e. excluding init time step for storage) time steps
    m_parameters['tm'] = m_parameters['timesteps'][1:]
    # within=m.t,
    # doc='Set of modelled timesteps'

    # support timeframes (e.g. 2020, 2030...)
    indexlist = set()
    for key in m_parameters['commodity_dict']["price"]:
        indexlist.add(tuple(key)[0])
    m_parameters['stf'] = pd.Index(indexlist, name='stf')
    # doc='Set of modeled support timeframes (e.g. years)'

    # site (e.g. north, middle, south...)
    indexlist = set()
    for key in m_parameters['commodity_dict']["price"]:
        indexlist.add(tuple(key)[1])
    m_parameters['sit'] = pd.Index(indexlist, name='sit')
    # doc='Set of sites'

    # commodity (e.g. solar, wind, coal...)
    indexlist = set()
    for key in m_parameters['commodity_dict']["price"]:
        indexlist.add(tuple(key)[2])
    m_parameters['com'] = pd.Index(indexlist, name='com')
    # doc='Set of commodities'

    # commodity type (i.e. SupIm, Demand, Stock, Env)
    indexlist = set()
    for key in m_parameters['commodity_dict']["price"]:
        indexlist.add(tuple(key)[3])
    m_parameters['com_type'] = pd.Index(indexlist, name='com_type')
    # doc='Set of commodity types'

    # process (e.g. Wind turbine, Gas plant, Photovoltaics...)
    indexlist = set()
    for key in m_parameters['process_dict']["inv-cost"]:
        indexlist.add(tuple(key)[2])
    m_parameters['pro'] = pd.Index(indexlist, name='pro')
    # doc='Set of conversion processes'

    # cost_type
    m_parameters['cost_type'] = pd.Index(m_parameters['cost_type_list'], name="cost_type")
    # doc='Set of cost types (hard-coded)'

    # tuple sets
    m_parameters['sit_tuples'] = pd.Index(m_parameters['site_dict']["area"].keys(),
                                          name="sit_tuples",
                                          tupleize_cols=False)
    # within=m.stf * m.sit
    # doc='Combinations of support timeframes and sites'
    m_parameters['com_tuples'] = pd.Index(m_parameters['commodity_dict']["price"].keys(),
                                          name="com_tuples",
                                          tupleize_cols=False)
    # within=m.stf * m.sit * m.com * m.com_type,
    # doc='Combinations of defined commodities, e.g. (2018,Mid,Elec,Demand)'
    m_parameters['pro_tuples'] = pd.Index(m_parameters['process_dict']["inv-cost"].keys(),
                                          name="pro_tuples",
                                          tupleize_cols=False)
    # within=m.stf * m.sit * m.pro,
    # doc='Combinations of possible processes, e.g. (2018,North,Coal plant)'
    m_parameters['com_stock'] = commodity_subset(m_parameters['com_tuples'], 'Stock')
    # within=m.com,
    # doc='Commodities that can be purchased at some site(s)')

    if m_parameters['mode']['int']:
        # tuples for operational status of technologies
        m_parameters['operational_pro_tuples'] = pd.Index(
                [(sit, pro, stf, stf_later) for (sit, pro, stf, stf_later)
                 in op_pro_tuples(m_parameters['pro_tuples'], m_parameters)],
                name="operational_pro_tuples", tupleize_cols=False)
        # within=m.sit * m.pro * m.stf * m.stf,
        # doc='Processes that are still operational through stf_later'
        #     '(and the relevant years following), if built in stf'
        #     'in stf.')

        # tuples for rest lifetime of installed capacities of technologies
        m_parameters['inst_pro_tuples'] = pd.Index([(sit, pro, stf)
                             for (sit, pro, stf)
                             in inst_pro_tuples(m_parameters)],
                             name="inst_pro_tuples", tupleize_cols=False)
        # within=m.sit * m.pro * m.stf,
        # doc='Installed processes that are still operational through stf')

    # commodity type subsets
    m_parameters['com_supim'] = commodity_subset(m_parameters['com_tuples'], 'SupIm')
    # within=m.com,
    # doc='Commodities that have intermittent (timeseries) input')
    m_parameters['com_demand'] = commodity_subset(m_parameters['com_tuples'], 'Demand')
    # within=m.com,
    # doc='Commodities that have a demand (implies timeseries)')
    m_parameters['com_env'] = commodity_subset(m_parameters['com_tuples'], 'Env')
    # within=m.com,
    # doc='Commodities that (might) have a maximum creation limit')

    # process tuples for area rule
    m_parameters['pro_area_tuples'] = pd.Index(
            m_parameters['proc_area_dict'].keys(), name="pro_area_tuples",
            tupleize_cols=False)
    # within=m.stf * m.sit * m.pro,
    # doc='Processes and Sites with area Restriction')

    # process input/output
    m_parameters['pro_input_tuples'] = pd.Index(
            [(stf, site, process, commodity)
             for (stf, site, process) in m_parameters['pro_tuples']
             for (s, pro, commodity) in tuple(m_parameters['r_in_dict'].keys())
             if process == pro and s == stf],
            name="pro_input_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro * m.com,
    # doc='Commodities consumed by process by site, e.g. (2020,Mid,PV,Solar)')
    m_parameters['pro_output_tuples'] = pd.Index(
            [(stf, site, process, commodity)
             for (stf, site, process) in m_parameters['pro_tuples']
             for (s, pro, commodity) in tuple(m_parameters['r_out_dict'].keys())
             if process == pro and s == stf],
            name="pro_output_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro * m.com,
    # doc='Commodities produced by process by site, e.g. (2020,Mid,PV,Elec)')

    # process tuples for maximum gradient feature
    m_parameters['pro_maxgrad_tuples'] = pd.Index([
        (stf, sit, pro) for (stf, sit, pro) in m_parameters['pro_tuples']
        if m_parameters['process_dict']['max-grad'][stf, sit, pro] < 1.0 / dt ],
        name="pro_maxgrad_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro,
    # doc='Processes with maximum gradient smaller than timestep length')

    # process tuples for partial feature
    m_parameters['pro_partial_tuples'] = pd.Index([
        (stf, site, process) for (stf, site, process) in m_parameters['pro_tuples']
        for (s, pro, _) in tuple(m_parameters['r_in_min_fraction_dict'].keys())
        if process == pro and s == stf ],
        name="pro_partial_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro,
    # doc='Processes with partial input')

    m_parameters['pro_partial_input_tuples'] = pd.Index([
        (stf, site, process, commodity)
        for (stf, site, process) in m_parameters['pro_partial_tuples']
        for (s, pro, commodity) in tuple(m_parameters['r_in_min_fraction_dict'].keys())
        if process == pro and s == stf ],
        name="pro_partial_input_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro * m.com,
    # doc='Commodities with partial input ratio, e.g. (2020,Mid,Coal PP,Coal)')

    m_parameters['pro_partial_output_tuples'] = pd.Index([
        (stf, site, process, commodity)
        for (stf, site, process) in m_parameters['pro_partial_tuples']
        for (s, pro, commodity) in tuple(m_parameters['r_out_min_fraction_dict'].keys())
        if process == pro and s == stf ],
        name="pro_partial_output_tuples", tupleize_cols=False)
    # within=m.stf * m.sit * m.pro * m.com,
    # doc='Commodities with partial input ratio, e.g. (Mid,Coal PP,CO2)')

    # Variables

    # A variable that equals 1. 
    # This is a hack to have linopy.LinearExpressions with constants in version
    # 0.0.8
    m.add_variables(coords=[], name="1")

    # costs
    m.add_variables(coords=[m_parameters['cost_type']], name="costs")
    # doc='Costs by type (EUR/a)')

    # commodity
    m.add_variables(lower=0,
                    coords=[m_parameters['tm'], m_parameters['com_tuples']],
                    name="e_co_stock")
    # within=pyomo.NonNegativeReals,
    # doc='Use of stock commodity source (MW) per timestep')

    # process
    m.add_variables(lower=0,
                    coords=[m_parameters['pro_tuples']],
                    name="cap_pro_new")
    # within=pyomo.NonNegativeReals,
    # doc='New process capacity (MW)')

    # process capacity as expression object
    # (variable if expansion is possible, else static)
    # TODO : To improve performance, replace the loop over index tuples
    # with operations on arrays
    # m.cap_pro = pyomo.Expression(
    #     m.pro_tuples,
    #     rule=def_process_capacity_rule,
    #     doc='Total process capacity (MW)')

    m_parameters['cap_pro'] = {
            (stf, sit, pro):def_process_capacity_rule(m,m_parameters,stf,sit,pro) 
            for (stf,sit,pro) in m_parameters['pro_tuples']}

    m.add_variables(lower=0,
                    coords=[m_parameters['t'], m_parameters['pro_tuples']],
                    name="tau_pro")
    # within=pyomo.NonNegativeReals,
    # doc='Power flow (MW) through process')

    m.add_variables(lower=0,
                    coords=[m_parameters['tm'],
                            m_parameters['pro_input_tuples']],
                    name="e_pro_in")
    # within=pyomo.NonNegativeReals,
    # doc='Power flow of commodity into process (MW) per timestep')

    m.add_variables(lower=0,
                    coords=[m_parameters['tm'],
                            m_parameters['pro_output_tuples']],
                    name="e_pro_out")
    # within=pyomo.NonNegativeReals,
    # doc='Power flow out of process (MW) per timestep')

    # Add additional features
    # called features are declared in distinct files in features folder
    if m_parameters['mode']['tra']:
        if m_parameters['mode']['dpf']:
            add_transmission_dc(m, m_parameters)
        else:
            add_transmission(m, m_parameters)
    if m_parameters['mode']['sto']:
        add_storage(m, m_parameters)
    if m_parameters['mode']['dsm']:
        add_dsm(m, m_parameters)
    if m_parameters['mode']['bsp']:
        add_buy_sell_price(m, m_parameters)
    if m_parameters['mode']['tve']:
        add_time_variable_efficiency(m, m_parameters)
    else:
        m_parameters['pro_timevar_output_tuples'] = set()
        # m.pro_timevar_output_tuples = pyomo.Set(
        #     within=m.stf * m.sit * m.pro * m.com,
        #     doc='empty set needed for (partial) process output')

    # Equation declarations
    # equation bodies are defined in separate functions, referred to here by
    # their name in the "rule" keyword.

    # Logic
    # This variable is set to 1 to allow constants in linear expressions
    m.add_constraint(m.variables["1"], "=", 1, name="const_1")

    # commodity
    m.res_vertex = pyomo.Constraint(
        m.tm, m.com_tuples,
        rule=res_vertex_rule,
        doc='storage + transmission + process + source + buy - sell == demand')
    m.res_stock_step = pyomo.Constraint(
        m.tm, m.com_tuples,
        rule=res_stock_step_rule,
        doc='stock commodity input per step <= commodity.maxperstep')
    m.res_stock_total = pyomo.Constraint(
        m.com_tuples,
        rule=res_stock_total_rule,
        doc='total stock commodity input <= commodity.max')
    m.res_env_step = pyomo.Constraint(
        m.tm, m.com_tuples,
        rule=res_env_step_rule,
        doc='environmental output per step <= commodity.maxperstep')
    m.res_env_total = pyomo.Constraint(
        m.com_tuples,
        rule=res_env_total_rule,
        doc='total environmental commodity output <= commodity.max')

    # process
    m.def_process_input = pyomo.Constraint(
        m.tm, m.pro_input_tuples - m.pro_partial_input_tuples,
        rule=def_process_input_rule,
        doc='process input = process throughput * input ratio')
    m.def_process_output = pyomo.Constraint(
        m.tm, (m.pro_output_tuples - m.pro_partial_output_tuples -
               m.pro_timevar_output_tuples),
        rule=def_process_output_rule,
        doc='process output = process throughput * output ratio')
    m.def_intermittent_supply = pyomo.Constraint(
        m.tm, m.pro_input_tuples,
        rule=def_intermittent_supply_rule,
        doc='process output = process capacity * supim timeseries')
    m.res_process_throughput_by_capacity = pyomo.Constraint(
        m.tm, m.pro_tuples,
        rule=res_process_throughput_by_capacity_rule,
        doc='process throughput <= total process capacity')
    m.res_process_maxgrad_lower = pyomo.Constraint(
        m.tm, m.pro_maxgrad_tuples,
        rule=res_process_maxgrad_lower_rule,
        doc='throughput may not decrease faster than maximal gradient')
    m.res_process_maxgrad_upper = pyomo.Constraint(
        m.tm, m.pro_maxgrad_tuples,
        rule=res_process_maxgrad_upper_rule,
        doc='throughput may not increase faster than maximal gradient')
    m.res_process_capacity = pyomo.Constraint(
        m.pro_tuples,
        rule=res_process_capacity_rule,
        doc='process.cap-lo <= total process capacity <= process.cap-up')

    m.res_area = pyomo.Constraint(
        m.sit_tuples,
        rule=res_area_rule,
        doc='used process area <= total process area')

    m.res_throughput_by_capacity_min = pyomo.Constraint(
        m.tm, m.pro_partial_tuples,
        rule=res_throughput_by_capacity_min_rule,
        doc='cap_pro * min-fraction <= tau_pro')
    m.def_partial_process_input = pyomo.Constraint(
        m.tm, m.pro_partial_input_tuples,
        rule=def_partial_process_input_rule,
        doc='e_pro_in = '
            ' cap_pro * min_fraction * (r - R) / (1 - min_fraction)'
            ' + tau_pro * (R - min_fraction * r) / (1 - min_fraction)')
    m.def_partial_process_output = pyomo.Constraint(
        m.tm,
        (m.pro_partial_output_tuples -
            (m.pro_partial_output_tuples & m.pro_timevar_output_tuples)),
        rule=def_partial_process_output_rule,
        doc='e_pro_out = '
            ' cap_pro * min_fraction * (r - R) / (1 - min_fraction)'
            ' + tau_pro * (R - min_fraction * r) / (1 - min_fraction)')

    #if m.mode['int']:
    #    m.res_global_co2_limit = pyomo.Constraint(
    #        m.stf,
    #        rule=res_global_co2_limit_rule,
    #        doc='total co2 commodity output <= global.prop CO2 limit')

    # costs
    m.def_costs = pyomo.Constraint(
        m.cost_type,
        rule=def_costs_rule,
        doc='main cost function by cost type')

    # objective and global constraints
    if m.obj.value == 'cost':
        m.res_global_co2_limit = pyomo.Constraint(
            m.stf,
            rule=res_global_co2_limit_rule,
            doc='total co2 commodity output <= Global CO2 limit')

        if m.mode['int']:
            m.res_global_co2_budget = pyomo.Constraint(
                rule=res_global_co2_budget_rule,
                doc='total co2 commodity output <= global.prop CO2 budget')

            m.res_global_cost_limit = pyomo.Constraint(
                m.stf,
                rule=res_global_cost_limit_rule,
                doc='total costs <= Global cost limit')

        m.objective_function = pyomo.Objective(
            rule=cost_rule,
            sense=pyomo.minimize,
            doc='minimize(cost = sum of all cost types)')

    elif m.obj.value == 'CO2':

        m.res_global_cost_limit = pyomo.Constraint(
            m.stf,
            rule=res_global_cost_limit_rule,
            doc='total costs <= Global cost limit')

        if m.mode['int']:
            m.res_global_cost_budget = pyomo.Constraint(
                rule=res_global_cost_budget_rule,
                doc='total costs <= global.prop Cost budget')
            m.res_global_co2_limit = pyomo.Constraint(
                m.stf,
                rule=res_global_co2_limit_rule,
                doc='total co2 commodity output <= Global CO2 limit')

        m.objective_function = pyomo.Objective(
            rule=co2_rule,
            sense=pyomo.minimize,
            doc='minimize total CO2 emissions')

    else:
        raise NotImplementedError("Non-implemented objective quantity. Set "
                                  "either 'cost' or 'CO2' as the objective in "
                                  "runme.py!")

    if dual:
        m.dual = pyomo.Suffix(direction=pyomo.Suffix.IMPORT)

    return m


# Constraints

# commodity

# vertex equation: calculate balance for given commodity and site;
# contains implicit constraints for process activity, import/export and
# storage activity (calculated by function commodity_balance);
# contains implicit constraint for stock commodity source term
def res_vertex_rule(m, tm, stf, sit, com, com_type):
    # environmental or supim commodities don't have this constraint (yet)
    if com in m.com_env:
        return pyomo.Constraint.Skip
    if com in m.com_supim:
        return pyomo.Constraint.Skip

    # helper function commodity_balance calculates balance from input to
    # and output from processes, storage and transmission.
    # if power_surplus > 0: production/storage/imports create net positive
    #                       amount of commodity com
    # if power_surplus < 0: production/storage/exports consume a net
    #                       amount of the commodity com
    power_surplus = - commodity_balance(m, tm, stf, sit, com)

    # if com is a stock commodity, the commodity source term e_co_stock
    # can supply a possibly negative power_surplus
    if com in m.com_stock:
        power_surplus += m.e_co_stock[tm, stf, sit, com, com_type]

    # if Buy and sell prices are enabled
    if m.mode['bsp']:
        power_surplus += bsp_surplus(m, tm, stf, sit, com, com_type)

    # if com is a demand commodity, the power_surplus is reduced by the
    # demand value; no scaling by m.dt or m.weight is needed here, as this
    # constraint is about power (MW), not energy (MWh)
    if com in m.com_demand:
        try:
            power_surplus -= m.demand_dict[(sit, com)][(stf, tm)]
        except KeyError:
            pass

    if m.mode['dsm']:
        power_surplus += dsm_surplus(m, tm, stf, sit, com)

    return power_surplus == 0

# stock commodity purchase == commodity consumption, according to
# commodity_balance of current (time step, site, commodity);
# limit stock commodity use per time step


def res_stock_step_rule(m, tm, stf, sit, com, com_type):
    if com not in m.com_stock:
        return pyomo.Constraint.Skip
    else:
        return (m.e_co_stock[tm, stf, sit, com, com_type] <=
                m.dt * m.commodity_dict['maxperhour']
                [(stf, sit, com, com_type)])


# limit stock commodity use in total (scaled to annual consumption, thanks
# to m.weight)
def res_stock_total_rule(m, stf, sit, com, com_type):
    if com not in m.com_stock:
        return pyomo.Constraint.Skip
    else:
        # calculate total consumption of commodity com
        total_consumption = 0
        for tm in m.tm:
            total_consumption += (
                m.e_co_stock[tm, stf, sit, com, com_type])
        total_consumption *= m.weight
        return (total_consumption <=
                m.commodity_dict['max'][(stf, sit, com, com_type)])


# environmental commodity creation == - commodity_balance of that commodity
# used for modelling emissions (e.g. CO2) or other end-of-pipe results of
# any process activity;
# limit environmental commodity output per time step
def res_env_step_rule(m, tm, stf, sit, com, com_type):
    if com not in m.com_env:
        return pyomo.Constraint.Skip
    else:
        environmental_output = - commodity_balance(m, tm, stf, sit, com)
        return (environmental_output <=
                m.dt * m.commodity_dict['maxperhour']
                [(stf, sit, com, com_type)])


# limit environmental commodity output in total (scaled to annual
# emissions, thanks to m.weight)
def res_env_total_rule(m, stf, sit, com, com_type):
    if com not in m.com_env:
        return pyomo.Constraint.Skip
    else:
        # calculate total creation of environmental commodity com
        env_output_sum = 0
        for tm in m.tm:
            env_output_sum += (- commodity_balance(m, tm, stf, sit, com))
        env_output_sum *= m.weight
        return (env_output_sum <=
                m.commodity_dict['max'][(stf, sit, com, com_type)])


# process

# process capacity (for m.cap_pro Expression)
def def_process_capacity_rule(m, m_parameters, stf, sit, pro):
    if m_parameters['mode']['int']:
        if (sit, pro, stf) in m_parameters['inst_pro_tuples']:
            if (sit, pro, min(m_parameters['stf']))\
                    in m_parameters['pro_const_cap_dict']:
                cap_pro = (m_parameters['process_dict']
                                       ['inst-cap']
                                       [stf, sit, pro]
                           * m.variables["1"])
            else:
                cap_pro = (sum(m.variables['cap_pro_new']
                                .loc[[(stf_built, sit, pro)]]
                               for stf_built in m.stf
                               if (sit, pro, stf_built, stf)
                               in m_parameters['operational_pro_tuples'])
                            + m_parameters['process_dict']
                                      ['inst-cap']
                                      [(min(m_parameters['stf']), sit, pro)]
                            * m.variables["1"])
        else:
            cap_pro = sum(m.variables['cap_pro_new'].loc[[(stf_built, sit, pro)]]
                for stf_built in m_parameters['stf']
                if (sit, pro, stf_built, stf)
                in m_parameters['operational_pro_tuples'])
    else:
        if (sit, pro, stf) in m_parameters['pro_const_cap_dict']:
            cap_pro = (
                    m_parameters['process_dict']['inst-cap'][[(stf, sit, pro)]]
                    * m.variables["1"])
        else:
            cap_pro = (m.variables['cap_pro_new'].loc[[(stf, sit, pro)]] +
                       + m_parameters['process_dict']
                                     ['inst-cap']
                                     [(stf, sit, pro)]
                       * m.variables["1"])
    return cap_pro

# process input power == process throughput * input ratio


def def_process_input_rule(m, tm, stf, sit, pro, com):
    return (m.e_pro_in[tm, stf, sit, pro, com] ==
            m.tau_pro[tm, stf, sit, pro] * m.r_in_dict[(stf, pro, com)])


# process output power = process throughput * output ratio
def def_process_output_rule(m, tm, stf, sit, pro, com):
    return (m.e_pro_out[tm, stf, sit, pro, com] ==
            m.tau_pro[tm, stf, sit, pro] * m.r_out_dict[(stf, pro, com)])


# process input (for supim commodity) = process capacity * timeseries
def def_intermittent_supply_rule(m, tm, stf, sit, pro, coin):
    if coin in m.com_supim:
        return (m.e_pro_in[tm, stf, sit, pro, coin] ==
                m.cap_pro[stf, sit, pro] * m.supim_dict[(sit, coin)]
                [(stf, tm)] * m.dt)
    else:
        return pyomo.Constraint.Skip


# process throughput <= process capacity
def res_process_throughput_by_capacity_rule(m, tm, stf, sit, pro):
    return (m.tau_pro[tm, stf, sit, pro] <= m.dt * m.cap_pro[stf, sit, pro])


def res_process_maxgrad_lower_rule(m, t, stf, sit, pro):
    return (m.tau_pro[t - 1, stf, sit, pro] -
            m.cap_pro[stf, sit, pro] *
            m.process_dict['max-grad'][(stf, sit, pro)] * m.dt <=
            m.tau_pro[t, stf, sit, pro])


def res_process_maxgrad_upper_rule(m, t, stf, sit, pro):
    return (m.tau_pro[t - 1, stf, sit, pro] +
            m.cap_pro[stf, sit, pro] *
            m.process_dict['max-grad'][(stf, sit, pro)] * m.dt >=
            m.tau_pro[t, stf, sit, pro])


def res_throughput_by_capacity_min_rule(m, tm, stf, sit, pro):
    return (m.tau_pro[tm, stf, sit, pro] >=
            m.cap_pro[stf, sit, pro] *
            m.process_dict['min-fraction'][(stf, sit, pro)] * m.dt)


def def_partial_process_input_rule(m, tm, stf, sit, pro, coin):
    # input ratio at maximum operation point
    R = m.r_in_dict[(stf, pro, coin)]
    # input ratio at lowest operation point
    r = m.r_in_min_fraction_dict[stf, pro, coin]
    min_fraction = m.process_dict['min-fraction'][(stf, sit, pro)]

    online_factor = min_fraction * (r - R) / (1 - min_fraction)
    throughput_factor = (R - min_fraction * r) / (1 - min_fraction)

    return (m.e_pro_in[tm, stf, sit, pro, coin] ==
            m.dt * m.cap_pro[stf, sit, pro] * online_factor +
            m.tau_pro[tm, stf, sit, pro] * throughput_factor)


def def_partial_process_output_rule(m, tm, stf, sit, pro, coo):
    # input ratio at maximum operation point
    R = m.r_out_dict[stf, pro, coo]
    # input ratio at lowest operation point
    r = m.r_out_min_fraction_dict[stf, pro, coo]
    min_fraction = m.process_dict['min-fraction'][(stf, sit, pro)]

    online_factor = min_fraction * (r - R) / (1 - min_fraction)
    throughput_factor = (R - min_fraction * r) / (1 - min_fraction)

    return (m.e_pro_out[tm, stf, sit, pro, coo] ==
            m.dt * m.cap_pro[stf, sit, pro] * online_factor +
            m.tau_pro[tm, stf, sit, pro] * throughput_factor)


# lower bound <= process capacity <= upper bound
def res_process_capacity_rule(m, stf, sit, pro):
    return (m.process_dict['cap-lo'][stf, sit, pro],
            m.cap_pro[stf, sit, pro],
            m.process_dict['cap-up'][stf, sit, pro])


# used process area <= maximal process area
def res_area_rule(m, stf, sit):
    if m.site_dict['area'][stf, sit] >= 0 and sum(
        m.process_dict['area-per-cap'][st, s, p]
        for (st, s, p) in m.pro_area_tuples
            if s == sit and st == stf) > 0:
        total_area = sum(m.cap_pro[st, s, p] *
                         m.process_dict['area-per-cap'][st, s, p]
                         for (st, s, p) in m.pro_area_tuples
                         if s == sit and st == stf)
        return total_area <= m.site_dict['area'][stf, sit]
    else:
        # Skip constraint, if area is not numeric
        return pyomo.Constraint.Skip


# total CO2 output <= Global CO2 limit
def res_global_co2_limit_rule(m, stf):
    if math.isinf(m.global_prop_dict['value'][stf, 'CO2 limit']):
        return pyomo.Constraint.Skip
    elif m.global_prop_dict['value'][stf, 'CO2 limit'] >= 0:
        co2_output_sum = 0
        for tm in m.tm:
            for sit in m.sit:
                # minus because negative commodity_balance represents creation
                # of that commodity.
                co2_output_sum += (- commodity_balance(m, tm,
                                                       stf, sit, 'CO2'))

        # scaling to annual output (cf. definition of m.weight)
        co2_output_sum *= m.weight
        return (co2_output_sum <= m.global_prop_dict['value']
                                                    [stf, 'CO2 limit'])
    else:
        return pyomo.Constraint.Skip


# CO2 output in entire period <= Global CO2 budget
def res_global_co2_budget_rule(m):
    if math.isinf(m.global_prop_dict['value'][min(m.stf_list), 'CO2 budget']):
        return pyomo.Constraint.Skip
    elif (m.global_prop_dict['value'][min(m.stf_list), 'CO2 budget']) >= 0:
        co2_output_sum = 0
        for stf in m.stf:
            for tm in m.tm:
                for sit in m.sit:
                    # minus because negative commodity_balance represents
                    # creation of that commodity.
                    co2_output_sum += (- commodity_balance
                                       (m, tm, stf, sit, 'CO2') *
                                       m.weight *
                                       stf_dist(stf, m))

        return (co2_output_sum <=
                m.global_prop_dict['value'][min(m.stf), 'CO2 budget'])
    else:
        return pyomo.Constraint.Skip


# total cost of one year <= Global cost limit
def res_global_cost_limit_rule(m, stf):
    if math.isinf(m.global_prop_dict["value"][stf, "Cost limit"]):
        return pyomo.Constraint.Skip
    elif m.global_prop_dict["value"][stf, "Cost limit"] >= 0:
        return(pyomo.summation(m.costs) <= m.global_prop_dict["value"]
               [stf, "Cost limit"])
    else:
        return pyomo.Constraint.Skip


# total cost in entire period <= Global cost budget
def res_global_cost_budget_rule(m):
    if math.isinf(m.global_prop_dict["value"][min(m.stf), "Cost budget"]):
        return pyomo.Constraint.Skip
    elif m.global_prop_dict["value"][min(m.stf), "Cost budget"] >= 0:
        return(pyomo.summation(m.costs) <= m.global_prop_dict["value"]
               [min(m.stf), "Cost budget"])
    else:
        return pyomo.Constraint.Skip


# Costs and emissions
def def_costs_rule(m, cost_type):
    #Calculate total costs by cost type.
    #Sums up process activity and capacity expansions
    #and sums them in the cost types that are specified in the set
    #m.cost_type. To change or add cost types, add/change entries
    #there and modify the if/elif cases in this function accordingly.
    #Cost types are
    #  - Investment costs for process power, storage power and
    #    storage capacity. They are multiplied by the investment
    #    factors. Rest values of units are subtracted.
    #  - Fixed costs for process power, storage power and storage
    #    capacity.
    #  - Variables costs for usage of processes, storage and transmission.
    #  - Fuel costs for stock commodity purchase.

    if cost_type == 'Invest':
        cost = \
            sum(m.cap_pro_new[p] *
                m.process_dict['inv-cost'][p] *
                m.process_dict['invcost-factor'][p]
                for p in m.pro_tuples)
        if m.mode['int']:
            cost -= \
                sum(m.cap_pro_new[p] *
                    m.process_dict['inv-cost'][p] *
                    m.process_dict['overpay-factor'][p]
                    for p in m.pro_tuples)
        if m.mode['tra']:
            # transmission_cost is defined in transmission.py
            cost += transmission_cost(m, cost_type)
        if m.mode['sto']:
            # storage_cost is defined in storage.py
            cost += storage_cost(m, cost_type)
        return m.costs[cost_type] == cost

    elif cost_type == 'Fixed':
        cost = \
            sum(m.cap_pro[p] * m.process_dict['fix-cost'][p] *
                m.process_dict['cost_factor'][p]
                for p in m.pro_tuples)
        if m.mode['tra']:
            cost += transmission_cost(m, cost_type)
        if m.mode['sto']:
            cost += storage_cost(m, cost_type)
        return m.costs[cost_type] == cost

    elif cost_type == 'Variable':
        cost = \
            sum(m.tau_pro[(tm,) + p] * m.weight *
                m.process_dict['var-cost'][p] *
                m.process_dict['cost_factor'][p]
                for tm in m.tm
                for p in m.pro_tuples)
        if m.mode['tra']:
            cost += transmission_cost(m, cost_type)
        if m.mode['sto']:
            cost += storage_cost(m, cost_type)
        return m.costs[cost_type] == cost

    elif cost_type == 'Fuel':
        return m.costs[cost_type] == sum(
            m.e_co_stock[(tm,) + c] * m.weight *
            m.commodity_dict['price'][c] *
            m.commodity_dict['cost_factor'][c]
            for tm in m.tm for c in m.com_tuples
            if c[2] in m.com_stock)

    elif cost_type == 'Environmental':
        return m.costs[cost_type] == sum(
            - commodity_balance(m, tm, stf, sit, com) * m.weight *
            m.commodity_dict['price'][(stf, sit, com, com_type)] *
            m.commodity_dict['cost_factor'][(stf, sit, com, com_type)]
            for tm in m.tm
            for stf, sit, com, com_type in m.com_tuples
            if com in m.com_env)

    # Revenue and Purchase costs defined in BuySellPrice.py
    elif cost_type == 'Revenue':
        return m.costs[cost_type] == revenue_costs(m)

    elif cost_type == 'Purchase':
        return m.costs[cost_type] == purchase_costs(m)

    else:
        raise NotImplementedError("Unknown cost type.")


def cost_rule(m):
    return pyomo.summation(m.costs)


# CO2 output in entire period <= Global CO2 budget
def co2_rule(m):
    co2_output_sum = 0
    for stf in m.stf:
        for tm in m.tm:
            for sit in m.sit:
                # minus because negative commodity_balance represents
                # creation of that commodity.
                if m.mode['int']:
                    co2_output_sum += (- commodity_balance(m, tm, stf, sit, 'CO2') *
                                       m.weight * stf_dist(stf, m))
                else:
                    co2_output_sum += (- commodity_balance(m, tm, stf, sit, 'CO2') *
                                       m.weight)

    return (co2_output_sum)
