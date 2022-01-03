import math
import random
import numpy as np

def constrainedLogGauss(center, width, bounds):
    new_val = 10.0 ** ( random.gauss(math.log10(center), width) )
    if new_val < bounds[0]: new_val = bounds[0]
    if new_val > bounds[1]: new_val = bounds[1]
    return new_val

def constrainedGauss(center, width, bounds):
    new_val = random.gauss(center, width)
    if new_val < bounds[0]: new_val = bounds[0]
    if new_val > bounds[1]: new_val = bounds[1]

    return new_val

def constrainedUniform(center, width, bounds):
    new_val = random.uniform(center - width / 2.0, center + width / 2.0)
    if new_val < bounds[0]: new_val = bounds[0]
    if new_val > bounds[1]: new_val = bounds[1]

    return new_val

def elSurrogateGauss(center, surrogate_width, bounds):
    surrogate_center = (center if center >= 1.0 else 2.0 - 1.0 / center )
    surrogate_bounds = [(bound if bound >= 1.0 else 2.0 - 1.0 / bound ) for bound in bounds]
    new_surrogate = constrainedGauss(surrogate_center, surrogate_width, surrogate_bounds)

    return (new_surrogate if new_surrogate > 1.0 else 1.0 / (2.0 - new_surrogate))

def elSurrogateUniform(center, surrogate_width, bounds):
    surrogate_center = (center if center >= 1.0 else 2.0 - 1.0 / center )
    surrogate_bounds = [(bound if bound >= 1.0 else 2.0 - 1.0 / bound ) for bound in bounds]
    new_surrogate = constrainedUniform(surrogate_center, surrogate_width, surrogate_bounds)

    return (new_surrogate if new_surrogate > 1.0 else 1.0 / (2.0 - new_surrogate))
    

def varyUnitVector(start_point, distance, fixed_values = {}):
    point_to_vary = [start_point[i] for i in range(len(start_point)) if i not in fixed_values.keys()]
    perturbation = [random.uniform(-1.0,1.0) for elem in point_to_vary]
    while sum([elem ** 2.0 for elem in perturbation]) <= 1.0:
        perturbation = [random.uniform(-1.0,1.0) for elem in point_to_vary]

    new_point = [point_to_vary[i] + perturbation[i] * distance for i in range(len(point_to_vary)) ]
    for index in fixed_values.keys():
        new_point.insert(index, fixed_values[index]) 
    new_size = np.sqrt(sum([elem ** 2.0 for elem in new_point]))
    
    new_point = [elem / new_size for elem in new_point]

    return new_point

def varyUpperHemisphereVector(start_point, distance):
    new_point = varyVector(start_point, distance)
    new_point = new_point [0:len(new_point) -1] + [abs(new_point[-1])]
    return new_point 

