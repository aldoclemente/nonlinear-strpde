#!/bin/bash

docker run --rm -d -p 8787:8787 -v $(pwd)/../../simulations:/home/user/simulations --name rstudio -e PASSWORD=password aldoclemente/fdapde-docker:rstudio


# connect to http://localhost:8787
# username: user
# password: password
