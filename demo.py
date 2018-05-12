'''
Demo of 2D calculation of rotational stability
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
    -0.5 * halfbeam,
    0.0,
    0.5 * halfbeam,
    1 * halfbeam])
sectionY = np.array([
    0.0,
    -0.5 * depth,
    -1.0 * depth,
    -0.5 * depth,
    0.0])

section = np.vstack([sectionX, sectionY])

def com(x, y):
    A = integrate.cumtrapz(x=x, y=-y)[-1]
    
    g_x = y
    inner = np.multiply(x, -g_x)
    com_x = integrate.cumtrapz(x=x, y=inner)[-1] / A

    y2 = np.square(y)
    com_y = integrate.cumtrapz(x=x, y=-y2)[-1] / (2 * A)
    com = np.array([com_x, com_y])
    com = np.round(com, 8)
    return com

sec_com = com(sectionX, sectionY)

waterX = np.linspace(-halfbeam, halfbeam, 100)
waterY = np.zeros(waterX.shape)

mass = 1406
supported_mass = mass / length
supported_liters = supported_mass / 1000

saved_sections = []
center_vectors = []

for boat_rotation in np.linspace(-pi/4, 0.5, 9):
    br = boat_rotation
    boat_R = np.matrix([[cos(br), -sin(br)], [sin(br), cos(br)]])
    rotated_com = boat_R * np.matrix(sec_com).T
    rotated_section = np.array(boat_R * np.matrix(section))
    
    water_level = 0.0

    # test
    for i in range(0, 5):
        break
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

    submerged_section = rotated_section.copy()
    rot_com = com(submerged_section[0,:], submerged_section[1,:])
    print('compare at %.2f' % (br,))
    print(rotated_com)
    print('vs')
    print(rot_com)
    submerged_section[1,:] = np.clip(submerged_section[1,:], None, 0)

    rot_cob = com(submerged_section[0,:], submerged_section[1,:])
    # print(rot_cob)

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
    1/0

for cv in center_vectors:
    plt.plot(*cv)
plt.title('center bouyancy vectors')
plt.show()
