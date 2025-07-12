from app import app

@app.route('/')
@app.route('/test')
def index():
    return "Hello, World!"