import numpy as np
import scipy as sp

from scipy.integrate import cumtrapz

# meters and kilograms

class Displacement2D(object):
    def __init__(self, area):
        self.area = area
    def as_water_mass(self, thickness=1):
        return 1000 * thickness * self.area

class Displacement3D(object):
    def __init__(self, volume):
        self.volume = volume
    def as_water_mass(self):
        return 1000 * self.volume

class Boat(object):
    X = 0
    Y = 1
    Z = 2
    def __init__(self):
        pass

    def to_points(self):
        raise NotImplementedError()

    def displacement(self):
        raise NotImplementedError()

class Boat2D(Boat):
    def displacement(self):
        vertices = self.to_points()
        sides = len(vertices)

        x_min_i, _ = np.argmin(vertices, axis=0)
        vertices = np.roll(vertices, -x_min_i, axis=0)
        x_max_i, _ = sides - np.argmax(np.flip(vertices), axis=0)

        bottoms = vertices[:x_max_i+1, :]
        tops = np.append(vertices[x_max_i-1:, :], vertices[0:1, :], axis=0)
        print('bottoms')
        print(bottoms)

        top_integral = abs(cumtrapz(tops[:, Boat.Y], tops[:, Boat.X],
            initial=0)[-1])
        bottom_integral = abs(cumtrapz(bottoms[:, Boat.Y], bottoms[:, Boat.X],
            initial=0)[-1])

        displacement_area = abs(top_integral - bottom_integral)
        print('displaced_area %.2f' % displacement_area)
        return Displacement2D(displacement_area)

        first_half = 0.0
        for i in range(0, n-2):
            first_half += vertices[i, X] * vertices[i + 1, Y] 

        middle = vertices[-1, X] * vertices[0, Y]

        second_half = 0.0
        for i in range(0, n-2):
            second_half += vertices[i + 1, X] * vertices[i, Y]

        end = vertices[1, X] * vertices[-1, Y]

        A = 0.5 * abs(first_half + middle - self_half - end)

class Boat3D(Boat):
    def displacement(self):
        return 0.0

class SimpleTrapezoidBoat(Boat2D):
    def to_points(self):
        x = [-2, -1, 1, 2]
        y = [0, -1, -1, 0]
        return np.array([x, y]).T

def will_it_float(boat, weight):
    return boat.displacement().as_water_mass() > weight

if __name__ == '__main__':
    print(will_it_float(SimpleTrapezoidBoat(), 90))
