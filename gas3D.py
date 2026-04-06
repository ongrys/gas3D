# Example file showing a circle moving on screen
import pygame, random, math
from time import time
import numpy as np

# length units are nanometers
# time units are nanoseconds
# velocity units are m/s
# mass unit is dalton internally and attogram when summed and displayed
# Energy is zepto J
# momentum is Dalton * m/s
# 1 liter is 1e-3 m3 ie. 1e24 nm3
# nano: -9, pico: -12, femto: -15, atto: -18, zepto: -21, yocto: -24


N = 500 
INIT_VEL = 500
NPART = 30

XS = 700                                # view window X pixels
YS = 700                                # view window Y pixels

NA = 6.0221415e23                       # Avogadro constant part/mole
DALTON = 1.66053906892e-27              # mass unit
Boltzmannk = 13.80649                   # Boltzmann constant x 10e24
vol_std = 0.0224                        # m3/mol
ppv = vol_std/NA                        # m3/particle
dx = ppv**(1/3)                         # m/particle
LEN = dx*(N**(1./3.))*1e9               # nanometers
                                        # modelled view is 0..LEN x 0..LEN x 0..LEN (z is parallel projection so brightness only)
pos = np.zeros((N,3),dtype=np.float64)
vel = np.zeros((N,3),dtype=np.float64)
mass = np.zeros((N),dtype=np.float64)
radius = np.zeros((N),dtype=np.float64)
color = np.zeros((N),dtype=int) # just hue
done = np.zeros((N),dtype=int)

dt = 0.001
czas = 0
fps = 0
totalU = 0
totalM = 0
maxV = 0
minV = 0

partitions = []

def assign_partitions():
    global partitions
    for i in range(NPART):
        for j in range(NPART):
            for k in range(NPART):
                partitions[i][j][k].clear()
    for i in range(N):
        try:
            partitions[int(pos[i][0]/LEN*NPART)][int(pos[i][1]/LEN*NPART)][int(pos[i][2]/LEN*NPART)].append(i)
        except IndexError:
            print(i,pos[i])
            raise KeyError

def draw_vel_spectrum(display_surface,maxV,num_buckets):
    buckets = []
    width = 500/num_buckets
    for i in range(num_buckets): buckets.append(0)
    for i in range(N):
        v = math.sqrt(np.dot(vel[i],vel[i]))
        j = int(v/maxV*num_buckets)
        if j > num_buckets - 1: j = num_buckets - 1
        buckets[j] += 1
    for i,v in enumerate(buckets):
        pygame.draw.rect(display_surface, "blue", pygame.Rect(750+i*width,700-v,width,v))


def display_stats(display_surface):
    font = pygame.font.Font('freesansbold.ttf', 12)

    text = font.render(f'Modelled dimension: {LEN:.2f} nm', True, "white", "black")
    display_surface.blit(text,(1080,10))
    text = font.render(f'Volume: {LEN**3/1e3:.5f} zeptoliters', True, "white", "black")
    display_surface.blit(text,(1080,25))
    text = font.render(f'Mass: {totalM:.5f} attograms', True, "white", "black")
    display_surface.blit(text,(1080,40))
    text = font.render(f'Internal energy: {totalU:.2f} zeptoJ', True, "white", "black")
    display_surface.blit(text,(1080,55))
    
    text = font.render(f'Particles: {N}', True, "white", "black")
    display_surface.blit(text,(920,10))
    text = font.render(f'Time: {czas:.1f} ns', True, "white", "black")
    display_surface.blit(text,(920,25))
    text = font.render(f'dt: {dt*1e3:.1f} ps', True, "white", "black")
    display_surface.blit(text,(920,40))
    text = font.render(f'FPS: {fps:.2f}', True, "white", "black")
    display_surface.blit(text,(920,55))

    text = font.render(f'Temperature: {totalU*1e3/N/Boltzmannk/3*2:.1f} K', True, "white", "black")
    display_surface.blit(text,(720,10))
    text = font.render(f'Energy density [m]: {totalU/totalM/1e3:.1f} kJ/kg', True, "white", "black") # zeptoJ /attog = J/kg
    display_surface.blit(text,(720,25))
    text = font.render(f'Energy density [V]: {totalU/LEN**3*1e3:.1f} J/liter', True, "white", "black") # zeptoJ /nm3 = J/ml
    display_surface.blit(text,(720,40))
    text = font.render(f'Density: {totalM/LEN**3*1e6:.5f} g/liter', True, "white", "black") # attog /nm3 --> kg/liter
    display_surface.blit(text,(720,55))

    text = font.render(f'Max V: {maxV:.1f} m/s', True, "white", "black")
    display_surface.blit(text,(720,70))
    text = font.render(f'Min V: {minV:.1f} m/s', True, "white", "black")
    display_surface.blit(text,(720,85))

    text = font.render(f'Press:', True, "white", "black")
    display_surface.blit(text,(1080,100))
    text = font.render(f'X for more Xenon', True, "white", "black")
    display_surface.blit(text,(1100,115))
    text = font.render(f'A for more Argon', True, "white", "black")
    display_surface.blit(text,(1100,130))

    draw_vel_spectrum(display_surface,2*INIT_VEL,20)
    text = font.render(f'{2*INIT_VEL:.1f} m/s', True, "white", "black")
    display_surface.blit(text,(1200,700))
    text = font.render("0", True, "white", "black")
    display_surface.blit(text,(750,700))
    
def calculate_stuff():
    global totalU,maxV,minV,totalM
    totalU = 0; totalM = 0
    minv2 = INIT_VEL/2; maxv2=INIT_VEL/2
    for i in range(N):
        V2 = vel[i][0]*vel[i][0]+vel[i][1]*vel[i][1]+vel[i][2]*vel[i][2]
        totalU += mass[i]*V2/2.0
        totalM += mass[i]
        if V2 > maxv2: maxv2 = V2
        if V2 < minv2: minv2 = V2
    maxV = math.sqrt(maxv2); minV = math.sqrt(minv2)
    totalU *= (DALTON * 1e21) # zeptoJ
    totalM *= (DALTON * 1e21) # atto gram

def add_particles(qty,m,r,h):
    global pos,vel,radius,mass,color,done,N,INIT_VEL
    bad_position = True
    while bad_position:
        new_pos = np.array([[random.random()*LEN*0.9+0.05*LEN,
                   random.random()*LEN*0.9+0.05*LEN,
                   random.random()*LEN*0.9+0.05*LEN]])
        bad_position = False
        for i in range(N):
            if math.sqrt((pos[i][0]-new_pos[0][0])**2
                        +(pos[i][1]-new_pos[0][1])**2
                        +(pos[i][2]-new_pos[0][2])**2) < radius[i]+r:
                bad_position = True
    new_vel = np.array([[random.random()-0.5,random.random()-0.5,random.random()-0.5]])
    len = math.sqrt(new_vel[0][0]**2 + new_vel[0][1]**2 + new_vel[0][2]**2)
    if len == 0: 
        new_vel[0] = 1.
        len = 1.
    new_vel = new_vel * (INIT_VEL/len)
    pos = np.append(pos,new_pos,axis=0)
    vel = np.append(vel,new_vel,axis=0)
    mass = np.append(mass,m)
    radius = np.append(radius,r)
    color = np.append(color,h)
    done = np.append(done,0)
    N+=1
        
    

def init_balls(vel0):
    # flat velocity profile
    for i in range(N):
        bad_position = True
        while bad_position:
            new_pos = [random.random()*LEN*0.9+0.05*LEN,
                       random.random()*LEN*0.9+0.05*LEN,
                       random.random()*LEN*0.9+0.05*LEN]
            bad_position = False
            for j in range(i):
                if math.sqrt((pos[j][0]-new_pos[0])**2
                            +(pos[j][1]-new_pos[1])**2
                            +(pos[j][2]-new_pos[2])**2) < 2*0.216:
                    bad_position = True
        pos[i] = new_pos
        new_vel = np.array([random.random()-0.5,random.random()-0.5,random.random()-0.5])
        len = math.sqrt(new_vel[0]**2 + new_vel[1]**2 + new_vel[2]**2)
        if len == 0: 
            new_vel[0] = 1.
            len = 1.
        new_vel = new_vel * (vel0/len)
        vel[i] = new_vel
        radius[i] = 0.031           # 216 pm for Xe, 31pm for Ar
        mass[i] = 39.8             # 131.3 for Xe, 39.8 for Ar
        color[i] = 240              # hue [0..360]
        
        
def screen_coords(actual_coords):
    return [10+XS/LEN*actual_coords[0],10+YS/LEN*actual_coords[1]]

def draw_radius(r):
    dr = XS/LEN*r
    if dr<3: return 3
    return dr

def draw_balls(screen):
    for i in range(N):
        clr = pygame.Color(0,0,0)
        L = int(pos[i][2]/LEN*100)
        clr.hsla = [color[i],L,L,20]
        pygame.draw.circle(screen, clr, screen_coords(pos[i]), draw_radius(radius[i]))

def calculate_collision_list(i,j,k):
    global partitions
    clist = []
    for x in [i-1,i,i+1]:
        for y in [j-1,j,j+1]:
            for z in [k-1,k,k+1]:
                if x < 0 or y < 0 or z < 0 or x == NPART or y == NPART or z == NPART:
                    continue
                clist.extend(partitions[x][y][z])
    return clist

def handle_collisions(i,clist,dt):
    global pos, vel, done
    for j in clist:
        if i == j: continue # no self-collisions
        if done[j]: continue # this has already collided
        dV = vel[i]-vel[j]
        dX = pos[i]-pos[j]
        A = np.dot(dV,dV)
        B = 2*np.dot(dV,dX)
        C = np.dot(dX,dX)-(radius[i]+radius[j])*(radius[i]+radius[j])
        DELTA = B*B-4*A*C
        if DELTA<=0: continue   # case of Delta == 0 is a glancing blow with no momentum transfer
        tc1 = (-B-math.sqrt(DELTA))/2./A
        tc2 = (-B+math.sqrt(DELTA))/2./A
        if tc1 < 0: 
            tc = tc2
        elif tc2 < 0: 
            tc = tc1
        else: 
            tc = min(tc1,tc2)
        if tc < 0: continue
        if tc <= dt:
            # print(f"Colliding in {tc} ns")
            # v1bar is the same as dV
            pbar = (dX+tc*dV)* (-1.0/(radius[i]+radius[j]))
            phi = mass[j]/mass[i]
            kc = 2./(phi+1)*np.dot(dV,pbar)
            # print(f"kc = {kc}")
            vic = dV-phi*kc*pbar+vel[j]     # velocities after the collision
            vjc = kc*pbar + vel[j]          # 
            # now move both to collision place
            pos[i] += vel[i]*tc
            pos[j] += vel[j]*tc
            # now finish the step
            vel[i] = vic;   vel[j] = vjc
            pos[i] += vel[i]*(dt-tc)
            pos[j] += vel[j]*(dt-tc)
            done[i] = 1;    done[j] = 1

            

def update_pos(dt):
    global pos, vel, done
    for i in range(NPART):
        for j in range(NPART):
            for k in range(NPART):
                for q in partitions[i][j][k]:
                    if done[q] == 0:
                        collision_list = calculate_collision_list(i,j,k)
                        handle_collisions(q,collision_list,dt)
                    if done[q] == 0:
                        pos[q] += vel[q]*dt
                        done[q] = 1
                    for d in range(3):
                        if pos[q][d] < 0: 
                            pos[q][d] *= -1.
                            vel[q][d] *= -1.
                        if pos[q][d] > LEN:
                            pos[q][d] = 2*LEN-pos[q][d]
                            vel[q][d] *= -1.
                            if 2*LEN-pos[q][d] < 0: print("Fuck")
                    # redundant check
                    for d in range(3):
                        if pos[q][d] < 0:   print("ZONX")
                        if pos[q][d] > LEN: print("ZONKY")
                            



# pygame setup
pygame.init()
screen = pygame.display.set_mode((1280, 720))
clock = pygame.time.Clock()
running = True
init_balls(INIT_VEL)
frame = 0
ts = time()
for i in range(NPART):
    partitions.append([])                      # x
    for j in range(NPART):
        partitions[i].append([])               # y
        for k in range(NPART):
            partitions[i][j].append([])        # z

while running:
    # poll for events
    # pygame.QUIT event means the user clicked X to close your window
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # fill the screen with a color to wipe away anything from last frame
    screen.fill("black")
    pygame.draw.rect(screen,"gray",pygame.Rect(10,10,700,700,width=1))
    done[:] = 0
    assign_partitions()
    update_pos(dt)
    draw_balls(screen)
    te = time()
    fps = 1/(te-ts)
    ts = te
    czas += dt
    calculate_stuff()
    dt = LEN/NPART/maxV/2   # in ns
    display_stats(screen)
    pygame.display.flip()
    keys = pygame.key.get_pressed()
    if keys[pygame.K_x]:
        add_particles(1,131.3,0.216,120)
        
    if keys[pygame.K_a]:
        add_particles(1,39.8,0.031,240)



pygame.quit()