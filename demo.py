'''
Demo of 2D calculation of center of bouyancy
'''

import numpy as np
from math import pi, cos, sin
from scipy import integrate
from matplotlib import pyplot as plt

length = 6.1
depth = 1.0 # from deck
beam = 2.71
halfbeam = beam * 0.5

sectionX = np.array([
    -1 * halfbeam,
    -0.2 * halfbeam,
    0.0,
    0.2 * halfbeam,
    1 * halfbeam])
sectionY = np.array([
    0.0,
    -0.75 * depth,
    -1.0 * depth,
    -0.75 * depth,
    0.0])

section = np.vstack([sectionX, sectionY])

# Center of mass is center of mass equation for full section (right now)
com_x = 0

x = section[0,:].flatten()
y = section[1,:].flatten()
A = integrate.cumtrapz(x=x, y=-y)[-1]
y2 = np.square(y)
com_y = integrate.cumtrapz(x=x, y=-y2)[-1] / (2 * A)

1/0

boat_rotation = pi/8
br = boat_rotation
water_level = 0.0

waterX = np.linspace(-halfbeam, halfbeam, 100)
waterY = np.zeros(waterX.shape)

mass = 1406
supported_mass = mass / length
supported_liters = supported_mass / 1000

saved_sections = []
center_vectors = []

for boat_rotation in np.linspace(-0.5, 0.5, 9):
    br = boat_rotation
    boat_R = np.matrix([[cos(br), -sin(br)], [sin(br), cos(br)]])
    rotated_section = np.array(boat_R * np.matrix(section))
    
    water_level = 0.0

    # test
    for i in range(0, 5):
        original_displace = integrate.cumtrapz(
            y=rotated_section[1].flatten(),
            x=rotated_section[0].flatten()) # ~= -2 m^3
        displacement = -original_displace[-1]

        if abs(displacement - supported_liters) < 0.0001:
            break
        delta = (displacement - supported_liters) / beam
        water_level += delta
        rotated_section[1,:] += delta

        # 1000 L per cubic meter
        # 1 L of water -> 1 kg -> 9.81 m/s -> 9.81 N
        # J24 -> 1406 kg, 3100 lb
        # 1 lbf ->  4.4482 N
    
    center_idx = len(rotated_section[0].flatten()) // 2


    fig = plt.figure()
    ax = fig.add_subplot('111')
    ax.set_title('Rotation: %.2f' % (br,))
    ax.plot(waterX, waterY)
    ax.plot(rotated_section[0,:], rotated_section[1,:])
    cv = (
        [0, rotated_section[0, center_idx]],
        [water_level, rotated_section[1, center_idx]],
    )
    center_vectors.append(cv)
    ax.plot(*cv,
        marker='x')
    ax.set_aspect(1.0)
    plt.tight_layout()
    plt.show()

for cv in center_vectors:
    plt.plot(*cv)
plt.title('center bouyancy vectors')
plt.show()
