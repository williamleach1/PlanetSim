import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astroquery.jplhorizons import Horizons
import datetime

body_ids = {
    "Mercury": 199,
    "Venus": 299,
    "Earth": 399,
    "Mars": 499,
    "Jupiter": 599,
    "Saturn": 699,
    "Uranus": 799,
    "Neptune": 899,
    "Pluto": 999,
    "Moon": 301,
    "Phobos": 401,
    "Deimos": 402,
    "Io": 501,
    "Europa": 502,
    "Ganymede": 503,
    "Callisto": 504,
    "Titan": 606
    # Add more moons if you want...
}


class Body:
    def __init__(self, mass, position, velocity):
        self.mass = mass
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)

solar_bodies = []

# # Gravitational constant in km^3 / kg / s^2
# G = 6.67430e-20

# for name, id in body_ids.items():
#     obj = Horizons(id=id, location='@ssb', epochs={'start':datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), 'stop':'2023-07-18', 'step':'1d'})
#     eph = obj.ephemerides()
#     #gm = obj.elements()['GM'][0]  # get GM
#     mass = 1#gm / G  # calculate mass
#     position = [eph['x'][0], eph['y'][0], eph['z'][0]]  # get position
#     velocity = [eph['vx'][0], eph['vy'][0], eph['vz'][0]]  # get velocity
#     solar_bodies.append(Body(mass, position, velocity))

class System:
    def __init__(self, bodies):
        self.bodies = bodies
        self.G = 1  # in our chosen units

    def compute_acceleration(self, body1, body2):
        r = body2.position - body1.position
        distance = np.linalg.norm(r)
        force = self.G * body1.mass * body2.mass / distance**3 * r
        acceleration = force / body1.mass
        return acceleration

    def update(self, dt):
        # Compute all forces
        forces = []
        for body1 in self.bodies:
            total_force = np.zeros(2)
            for body2 in self.bodies:
                if body1 != body2:
                    total_force += self.compute_acceleration(body1, body2)
            forces.append(total_force)

        # Update velocities
        for i in range(len(self.bodies)):
            self.bodies[i].velocity += forces[i] * dt

        # Update positions
        for body in self.bodies:
            body.position += body.velocity * dt
            print(body.position)

    def animate(self, anchor_body=None):
        fig, ax = plt.subplots()

        def init():
            return []

        def update(frame):
            ax.clear()
            self.update(0.0001)
            markers = []
            for body in self.bodies:
                marker, = ax.plot(*body.position, 'o')
                markers.append(marker)
            if anchor_body is not None:
                ax.set_xlim(-0.2, 0.5)
                ax.set_ylim(-0.2, 0.2)
            else:
                ax.set_xlim(-10, 10)
                ax.set_ylim(-10, 10)
            return markers

        ani = animation.FuncAnimation(fig, update, init_func=init, blit=True, interval=20, frames=1000)

        plt.show()

# Initial conditions
earth = Body(1, [0, 0], [0, -0.1])  # Earth is at rest at the origin
moon = Body(0.0123, [0.0, 3*0.0257], [0.5*6.5, 0.0])  # Moon has initial velocity to keep it in orbit
moon2 = Body(0.00123, [0.0, -3.3*0.0257], [0.5*-9, 0.0])
moon3 = Body(0.0123, [0.0, 4*0.0257], [0.1*9, 0.0])
moon4 = Body(0.05, [0.0, -4*0.0257], [-0.35*9, 0.0])
# Create system and animate§
system = System([earth, moon2, moon4])
system.animate(earth)  # using the first body (Earth in this case) as the anchor
