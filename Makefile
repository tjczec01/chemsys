init:
    pip install -r requirements.txt

test:
    
    pip install -r requirements_test.txt
    nosetests tests
