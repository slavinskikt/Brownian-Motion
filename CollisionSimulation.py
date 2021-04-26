#2D particle collision program
import numpy as np
import random as rand
from matplotlib import pyplot as plt
from matplotlib import animation

class Ball():
    def __init__(self, vx, vy, x, y, a, mass, radius, color):
        self.vx = vx
        self.vy = vy
        self.x = x
        self.y = y
        self.a = a
        self.mass = mass
        self.radius = radius
        self.color = color
        self.top = self.y + self.radius
        self.bottom = self.y - self.radius
        self.left = self.x - self.radius
        self.right = self.x + self.radius
        
    def move(self, dt):
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.vx += self.a*dt
        self.vy += self.a*dt
        
        self.top = self.y + self.radius
        self.bottom = self.y - self.radius
        self.left = self.x - self.radius
        self.right = self.x + self.radius
        
class Box():
    def __init__(self, length, height):
        self.length = length
        self.height = height
        self.top = height
        self.bottom = 0
        self.right = length
        self.left = 0
        
def simulateCollision(N = 50,
                    velocity = 2,
                    big_radius = 0.25,
                    small_radius = 0.05,
                    big_mass = 1,
                    small_mass = 0.9,
                    time = 5, 
                    FPS = 60, 
                    length = 5, 
                    height = 5):
    '''
    This program models a physical phenomena known as brownian motion. 
    Brownian motion is a process that governs diffusion processes. The
    kinetics of water molecules, for example, create a movement pattern that
    appears random. This phenomena was discovered in 1827 by Robert Brown after
    observing the random movement fluctuations of a pollen grain in water.
    This is what this program aims to simulate.
    
    To accomplish this, a physics engine was created. Specifically, collision
    dynamics were modeled in 2D, in order to achieve the desired effect. The
    collision dynamics were considered to be perfectly elastic.
    
    As arguments, this function (simulateCollision) takes the following:
        -N = the number of particles aside from the pollen grain
        -velocity = the maximum starting velocity of a given molecule
        -big_radius = the radius of the 'pollen grain'
        -small_radius = the radius of a 'water molecule'
        -big_mass = the mass of the 'pollen grain'
        -small_mass = the mass of the 'water molecule'
        -time = the approximate duration of the simulation
        -FPS = the desired frames per second of the simulation
        -length = the length of the viewing window
        -height = the height of the viewing window
    
    This program utilizes balls to model the aforementioned water and pollen.
    A Ball class was created to achieve this.
    '''
    
    def mover(dt, time, box, balls, Xdata, Ydata):
        #moves the balls in accordance with collision dynamics
        f = 0
        frames = 1/dt*time
        while f < frames:
            j = 0
            for ball in balls:
                handleBoxCollision(ball, box)
                k = 0
                for ball2 in balls:
                    if k != j:
                        handleBallCollision(ball, ball2)
                    k += 1
                        
                ball.move(dt)
    
                # print(ball.x, ball.y)
                Xdata[f, j] = ball.x 
                Ydata[f, j] = ball.y
                j += 1
            f += 1
        return Xdata, Ydata


    def handleBoxCollision(ball, box):
        #algorithm for checking collisions with box walls
        if ball.top >= box.top or ball.bottom <= box.bottom:
            ball.vy = -ball.vy
        if ball.right >= box.right or ball.left <= box.left:
            ball.vx = -ball.vx
        # return ball.vx, ball.vy
        
    def handleBallCollision(ball1, ball2):
        #algorithm for checking inter-ball collisions
        def calculateVelocity(v1i, v2i, X1, X2, m1, m2):
            #calculates the final velocity vectors of the colliding balls
            v1f = v1i-(2*m2/(m1+m2))* \
                np.dot(v1i-v2i, X1-X2)/(np.linalg.norm(X1-X2)**2)*(X1-X2)
            return v1f
        
        dist = np.sqrt((ball2.x - ball1.x)**2 + (ball2.y - ball1.y)**2)
        R = ball1.radius + ball2.radius
        if dist <= R:
            X1 = np.array([ball1.x, ball1.y])
            X2 = np.array([ball2.x, ball2.y])
            v1i = np.array([ball1.vx, ball1.vy])
            v2i = np.array([ball2.vx, ball2.vy])
            m1 = ball1.mass
            m2 = ball2.mass
            ball1.vx, ball1.vy = calculateVelocity(v1i, v2i, X1, X2, m1, m2)
            ball2.vx, ball2.vy = calculateVelocity(v2i, v1i, X2, X1, m2, m1)
            
    def init():
        #initalizes the animation 
        for n in range(len(balls)):
            patch = circles[n]
            patch.center = (Xdata[0,n], Ydata[0,n])
            ax.add_patch(patch)
        return circles
        
    def animate(i):
        #animation function
        for n in range(len(balls)):
            patch = circles[n]
            patch.center = (Xdata[i,n], Ydata[i,n])
        return circles
    
    
    #time step
    dt = 1/FPS
    #total number of frames in animation
    totalframes = time*FPS
    
    #acceleration is assumed to be zero
    a = 0
    #created instance of Ball object for larger ball object (pollen grain)
    bigBall = Ball(0, 0, height/2, length/2, 
                   a, big_mass, big_radius, 'red')
    #create instance of the Box object to be used as the viewing window
    box = Box(length, height)
    
    #create instances of N smaller ball objects, as well as circle patches
    #to be used in the animation (which represent the balls), and store these
    #ball and circle objects in their respective arrays
    balls = []
    circles = []
    for _ in range(N):
        xi = rand.uniform(small_radius, box.length - small_radius)
        yi = rand.uniform(small_radius, box.height - small_radius)
        vxi = rand.uniform(-velocity, velocity)
        vyi = rand.uniform(-velocity, velocity)
        ball = Ball(vxi, vyi, xi, yi, a, small_mass, small_radius, 'blue')
        balls.append(ball)
        circles.append(plt.Circle((ball.x,ball.y), ball.radius,
                                  fc = ball.color))
    balls.append(bigBall)
    circles.append(plt.Circle((bigBall.x,ball.y), bigBall.radius,
                              fc = bigBall.color))
    
    #initialize arrays to store values of the x and y positions throughout the
    #duration of the simulation
    Xdata = np.zeros((len(np.arange(0,time,dt)), N+1))
    Ydata = np.zeros((len(np.arange(0,time,dt)), N+1))
    
    #simluate movement
    Xdata, Ydata = mover(dt, time, box, balls, Xdata, Ydata)
    
    #create figure objects for animatoin putposes
    fig, ax = plt.subplots()
    ax = plt.axes(xlim=(0, box.length), ylim=(0, box.height))
    
    ax.set_xlabel("Horizontal Displacement")   
    ax.set_ylabel("Vertical Displacement")
    ax.set_title("Brownian Motion Simulation")
    
    #animate
    ani = animation.FuncAnimation(fig, animate, 
                                  init_func = init, 
                                  frames = totalframes,
                                  interval = dt*1000,
                                  blit = True,
                                  cache_frame_data = True)
    plt.show()
    return Xdata, Ydata       
        

X, Y = simulateCollision(N = 50,
            velocity = 2,
            big_radius = 0.25,
            small_radius = 0.05,
            big_mass = 1,
            small_mass = 0.9,
            time = 10,
            FPS = 120,
            length = 5,
            height = 5)   
     
        
        