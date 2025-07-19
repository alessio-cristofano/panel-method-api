
from flask import Flask, Response, jsonify
from flask_cors import CORS

from app.api import solver_bp


# from .config import Config


def create_app():
    app = Flask(__name__)
    CORS(app)
    # app.config.from_object(Config)
    # app.teardown_appcontext(close_kv)

    # init_db(app)

    app.register_blueprint(solver_bp)

    return app