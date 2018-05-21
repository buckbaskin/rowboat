'''
Demo of 2D calculation of rotational stability

Steps:

Build section
    Add hull properties
Visualize section
Calculate COM (should rotate with hull because all rigid)
build wave profile

rotate
find submerged level
find displaced COM
find lever arm

Visualize rotation profile

Quantify rotational stability compared to wind resistance
    Quantify sail area and lever arm
    Pick wind speed -> knocking force
    compare to righting force

Bless sailboat as stable or not yet
'''

import numpy as np
from math import pi, cos, sin
from scipy import integrate
from matplotlib import pyplot as plt

length = 6.1
depth = 1.0 # from deck
beam = 2.71
halfbeam = beam * 0.5
mass = 1406
freeboard = 0.9 # m
draft = 1.22 # m
keel_mass = 431 # kg
keel_depth = freeboard + draft # m
keel_com = -(depth + (keel_depth - depth) / 2) # (fin keel, no division for bulb keel)
hull_mass = mass - keel_mass
supported_mass = hull_mass / length
supported_liters = supported_mass / 1000

wind_velocity = 10 # mph
wind_velocity = wind_velocity * 0.447 # m / sec
wind_density = 1.225 # kg / m^3
sail_area_100p = 24.26 # m^2
mast_height = 10 # m (approx)
mast_lever = mast_height / 3

# wind_knock = 

sectionX = np.array([
    -1 * halfbeam,
    -0.95 * halfbeam,
    -0.8 * halfbeam,
    -0.15 * halfbeam,
    0.0,
    0.15 * halfbeam,
    0.8 * halfbeam,
    0.95 * halfbeam,
    1 * halfbeam])
sectionY = np.array([
    0.0,
    -0.4 * depth,
    -0.55 * depth,
    -0.85 * depth,
    -1.0 * depth,
    -0.85* depth,
    -0.55 * depth,
    -0.4 * depth,
    0.0])

visualize_section = False
if visualize_section:
    fig = plt.figure()
    ax = fig.add_subplot('111')
    ax.plot(sectionX, sectionY)
    ax.plot([0, 0, 0], [0, keel_com, -keel_depth], marker='x')
    ax.set_aspect(1.0)
    plt.show()
    import sys
    sys.exit(0)

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
keel_com = np.matrix([[0, keel_com]])
boat_com = (keel_com * keel_mass + sec_com * hull_mass) / (keel_mass + hull_mass)

waterX = np.linspace(-halfbeam, halfbeam, 100)
waterY = np.zeros(waterX.shape)

polarX = []
polarY = []
center_vectors = []

resolution = 101

max_deg = 30
max_rad = max_deg / 180 * pi

unstable_design = False

for boat_rotation in np.linspace(-max_rad, max_rad, resolution):
    br = boat_rotation
    boat_R = np.matrix([[cos(br), -sin(br)], [sin(br), cos(br)]])
    rotated_com = boat_R * np.matrix(boat_com).T
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
        rotated_com[1,:] += delta
        rotated_section[1,:] += delta

        # 1000 L per cubic meter
        # 1 L of water -> 1 kg -> 9.81 m/s -> 9.81 N
        # J24 -> 1406 kg, 3100 lb
        # 1 lbf ->  4.4482 N

    submerged_section = rotated_section.copy()
    submerged_section[1,:] = np.clip(submerged_section[1,:], None, 0)

    rot_cob = com(submerged_section[0,:], submerged_section[1,:])

    righting_lever = np.matrix(rot_cob).T - rotated_com
    polarX.append(br * (180 / pi)) # degrees
    righting_moment = float(righting_lever[0]) * (mass * 9.81) # Newton meters
    if abs(br) > 0.00001:
        righting_moment *= (-br / abs(br))
    righting_torque = righting_moment / 1.3556 # ft-lb
    polarY.append(righting_torque)

    rcomx, rcomy = rotated_com
    rcobx, rcoby = rot_cob

    if resolution <= 7 or (righting_moment < 0 and unstable_design == False):
        fig = plt.figure()
        ax = fig.add_subplot('111')
        title = 'Ship Display (%.1f deg)' % (br * 180 / pi,)
        if righting_moment < 0:
            title += ' Unstable'
        ax.set_title(title)
        ax.plot(waterX, waterY)
        ax.plot(rotated_section[0,:], rotated_section[1,:])
        cv = (
            [rcomx, rcobx],
            [rcomy, rcoby],
        )
        center_vectors.append(cv)
        ax.plot(*cv,
            marker='x')
        ax.set_aspect(1.0)
        plt.tight_layout()
        plt.show()
    if righting_moment <= 0:
        if br > 0 or br < 0:
           unstable_design = True

fig = plt.figure()
ax = fig.add_subplot('111')
ax.plot(polarX, polarY)
ax.set_title('Stability Curve')
plt.show()

print('Stable Design? %s' % (not unstable_design,))
