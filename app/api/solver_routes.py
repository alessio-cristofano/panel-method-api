from app import app

@app.route('/')
@app.route('/solve')
def index():
    return "Hello, World!"