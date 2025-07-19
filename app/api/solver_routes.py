from app.solver import solve

from flask import jsonify,Blueprint,request

solver_bp = Blueprint("solver", __name__)


@solver_bp.route("/solve", methods=["POST"])
def solve_model():
    solver_input = request.get_json()
    print(solver_input)
    dict_output = solve(solver_input)
    return(jsonify(dict_output),200)