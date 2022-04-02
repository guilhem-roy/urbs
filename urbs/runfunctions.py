import os
import pyomo.environ
from pyomo.opt.base import SolverFactory
import linopy
from datetime import datetime, date
from .model import create_model
from .report import *
from .plot import *
from .input import *
from .validation import *
from .saveload import *


def prepare_result_directory(result_name):
    """ create a time stamped directory within the result folder.

    Args:
        result_name: user specified result name

    Returns:
        a subfolder in the result folder 
    
    """
    # timestamp for result directory
    now = datetime.now().strftime('%Y%m%dT%H%M')

    # create result directory if not existent
    result_dir = os.path.join('result', '{}-{}'.format(result_name, now))
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    return result_dir


def setup_solver(solver, logfile='solver.log'):
    """ """
    options = dict()
    if solver == 'gurobi':
        # reference with list of option names
        # http://www.gurobi.com/documentation/5.6/reference-manual/parameters
        options["logfile"] = logfile
        # options["timelimit"] = 7200  # seconds
        # options["mipgap"] = 5e-4  # default = 1e-4
    elif solver == 'glpk':
        # reference with list of options
        # execute 'glpsol --help'
        options["log"]=logfile
        # options["tmlim"] = 7200  # seconds
        # options["mipgap"] = .0005
    elif solver == 'cplex':
        options["log"]=logfile
    else:
        print("Warning from setup_solver: no options set for solver "
              "'{}'!".format(optim.name))
    return options


def run_scenario(input_files, solver, timesteps, scenario, result_dir, dt,
                 objective, plot_tuples=None,  plot_sites_name=None,
                 plot_periods=None, report_tuples=None,
                 report_sites_name=None):
    """ run an urbs model for given input, time steps and scenario

    Args:
        - input_files: filenames of input Excel spreadsheets
        - Solver: the user specified solver
        - timesteps: a list of timesteps, e.g. range(0,8761)
        - scenario: a scenario function that modifies the input data dict
        - result_dir: directory name for result spreadsheet and plots
        - dt: length of each time step (unit: hours)
        - objective: objective function chosen (either "cost" or "CO2")
        - plot_tuples: (optional) list of plot tuples (c.f. urbs.result_figures)
        - plot_sites_name: (optional) dict of names for sites in plot_tuples
        - plot_periods: (optional) dict of plot periods
          (c.f. urbs.result_figures)
        - report_tuples: (optional) list of (sit, com) tuples
          (c.f. urbs.report)
        - report_sites_name: (optional) dict of names for sites in
          report_tuples

    Returns:
        the urbs model instance
    """

    # sets a modeled year for non-intertemporal problems
    # (necessary for consitency)
    year = date.today().year

    # scenario name, read and modify data for scenario
    sce = scenario.__name__
    data = read_input(input_files, year)
    data = scenario(data)
    validate_input(data)
    validate_dc_objective(data, objective)

    # create model
    (prob, prob_params) = create_model(data, dt, timesteps, objective)
    # prob_filename = os.path.join(result_dir, 'model.lp')
    # prob.write(prob_filename, io_options={'symbolic_solver_labels':True})

    # refresh time stamp string and create filename for logfile
    log_filename = os.path.join(result_dir, '{}.log').format(sce)

    # solve model and read results
    solver_options = setup_solver(solver, log_filename)
    prob.solve(solver, log_fn=log_filename, **solver_options)
    print(prob)

    # save problem solution (and input data) to HDF5 file
    save(prob, prob_params, os.path.join(result_dir, '{}.h5'.format(sce)))

    # write report to spreadsheet
    report(
        prob,
        prob_params,
        os.path.join(result_dir, '{}.xlsx').format(sce),
        report_tuples=report_tuples,
        report_sites_name=report_sites_name)

    # result plots
    result_figures(
        prob,
        prob_params,
        os.path.join(result_dir, '{}'.format(sce)),
        timesteps,
        plot_title_prefix=sce.replace('_', ' '),
        plot_tuples=plot_tuples,
        plot_sites_name=plot_sites_name,
        periods=plot_periods,
        figure_size=(24, 9))

    return prob
