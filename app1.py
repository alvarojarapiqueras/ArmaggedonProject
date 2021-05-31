from flask import Flask, render_template, url_for,session, request, \
    copy_current_request_context
from flask_socketio import SocketIO, emit, join_room, leave_room, \
    close_room, rooms, disconnect
import numpy as np
import json
 


app = Flask(__name__)
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

@app.route('/')
def index():
    return render_template('index.html', **locals())
   

@socketio.event
def my_event(message):

    x0= float(message['x0'])
    y0= float(message['y0'])
    z0= float(message['z0'])
    v0= float(message['v0'])

    class posInicial:
        def __init__(self, x0, y0, z0, v, n_frames, R_anim):
            self.Rt=6378    #km
            self.Rp = 30000 #r perigeo trayectoria misil, altura de impacto
            self.x0=x0
            self.z0=z0
            self.y0=y0
            self.v=v    #v meteo
            self.n_frames=n_frames #numero de frames para dar una vuelta a la Tierra
            self.R_anim=R_anim
            self.pos0=np.array([x0, y0, z0])

        def inicial(self):
            
            pos_unitario = self.pos0/np.linalg.norm(self.pos0)

            return pos_unitario

        def vect_pos(self):
            
            
            r0 = posInicial.distancia(self)
            pos_unitario = posInicial.inicial(self)
            vect = posInicial.tiempos_frames(self)
            tf = vect[0]
            R = self.R_anim/6378 # 1km realidad = 1/6378 en la animacion
            pos = np.zeros((tf,3))

            for i in range(tf):
                d = r0/(tf-1)
                pos[i]= self.pos0 - i*d*pos_unitario

            pos_impacto = pos_unitario

            pos_anim = pos * R #vector con posiciones en la animacion

            r=np.zeros((tf,1))
            for i in range(tf):
                r[i]=np.sqrt(pos_anim[i,0]**2+pos_anim[i,1]**2+pos_anim[i,2]**2)

            return pos_anim, r, pos_impacto

        def coord(self):

            if (self.z0 >= 0):
                n=1
                q=-1
            else:
                n=-1
                q=1
            if (self.y0 >=0):
                m=1
            else:
                m=-1

            if (self.x0 >=0):
                l=0
            else:
                l=np.pi/2
        
            lat = np.arcsin(np.absolute(self.pos0[2])/np.linalg.norm(self.pos0)) * n   #rad

            beta = np.pi/2 + np.arcsin(np.absolute(self.pos0[2])/np.linalg.norm(self.pos0))*q

            #longitud= angulo con meridiano greenwich, plano (0,1,0)

            long = (np.arcsin(np.absolute(self.y0)/(np.linalg.norm(self.pos0)*np.cos(lat))) + l) * m   #rad

            vect = [lat, beta, long]

            return vect

        def distancia(self):
            dist = np.linalg.norm(self.pos0)
            dist = dist - self.Rp

            return dist

        def tiempos_frames(self):

            h = self.n_frames/24  #1 hora = 100/24 frames
            t = posInicial.distancia(self)/self.v #tiempo en s
            th = t/3600 #h
            tf = th* h  # numero de frames en una simulacion
            tf = int(tf) #numero de frames = entero inmediatamente inferior
            vect = [tf, h]
            return vect
        



    ################################################################################
    ################################################################################

    class misil:
        def __init__(self,  R_anim, vect,r_meteo,lat_lon):

            # vect es h y tf;
            # r_meteo es dist meteo (r)
            # lat_lont es lat beta lon

            self.Rt= 6378    #km
            self.R_anim = R_anim
            self.v0= 10 #km/s
            self.muT = 3.986*10**5    #km^3/s^2
            self.r0 = self.Rt
            self.Rp = 30000 #km Radio del apogeo
            self.Q = self.r0*self.v0**2/self.muT #conocemos
            self.a = self.r0/(2-self.Q) #conocemos
            self.e = self.Rp/self.a-1
            self.phi0 = np.arccos(np.sqrt((1-self.e**2)/(self.Q*(2-self.Q))))
            self.tf = vect[0]
            self.h = vect[1]
            self.r_meteo = r_meteo
            self.lat = lat_lon[0]
            self.beta = lat_lon[1]
            self.long = lat_lon[2]

        def ang_alcance(self):
            
            denominador = np.sqrt(1-self.Q*(2-self.Q)*(np.cos(self.phi0))**2)
            numerador = (1-self.Q*(np.cos(self.phi0)**2))
            alcance=np.arccos(numerador/denominador)

            return alcance

        def t_impacto(self):

            E0=np.arccos((self.e-np.cos(misil.ang_alcance(self)))/(1-self.e*np.cos(misil.ang_alcance(self))))
            tv = 2*np.sqrt(self.a**3/self.muT)*(np.pi-E0+self.e+self.e*np.sin(E0)) #s
            tv_hora=tv/3600 #Tiempo en horas hasta impacto

            n_frames_hasta_impacto = tv_hora * self.h   #h es la relacion frames/hora
            n_frames_hasta_impacto = int(n_frames_hasta_impacto) #n frames en los que tenemos que calcular la trayectoria
                                                        #del misil
            
            

                
            n_frames_hasta_misil = self.tf - n_frames_hasta_impacto #n frames hasta que se lanza el misil        
            n_frames = [n_frames_hasta_impacto, n_frames_hasta_misil, self.tf]


            return n_frames
        
        def rmisil(self):

            pmisil = self.r0*self.Q*(np.cos(self.phi0))**2
            nu0 = np.pi - misil.ang_alcance(self)
            r_misil = np.zeros((misil.t_impacto(self)[0],1))
            r_misil = np.zeros((misil.t_impacto(self)[2],1))
            nu = np.linspace(nu0,np.pi,num=misil.t_impacto(self)[0])

        
            for i in range(misil.t_impacto(self)[2]):
                if (i < misil.t_impacto(self)[1]):
                    r_misil[i]=0
                    
                else:
                    r_misil[i]= pmisil/(1+self.e*np.cos(nu[i-misil.t_impacto(self)[1]]))

            
                    
            return r_misil, nu
        
        def coord_misil(self):
            coord_cart = np.zeros((misil.t_impacto(self)[2],3))
            alpha = np.zeros((misil.t_impacto(self)[0],1))

            r_misil , nu = misil.rmisil(self)
            n_frames = misil.t_impacto(self)

            for i in range(n_frames[2]):
                if (i < n_frames[1]):
                    pass
                else:
                    alpha[i-n_frames[1]] = self.long - misil.ang_alcance(self) + (nu[i-n_frames[1]] - nu[0])
                    coord_cart[i,0] = r_misil[i] * np.sin(self.beta) * np.cos(alpha[i-n_frames[1]])
                    coord_cart[i,1] = r_misil[i] * np.sin(self.beta) * np.sin(alpha[i-n_frames[1]])
                    coord_cart[i,2] = r_misil[i] * np.cos(self.beta)

            R = self.R_anim/6378 # 1km realidad = 1/6378 en la animacion
            cart_anim = coord_cart * R #vector con posiciones en la animacion
            

            return cart_anim, coord_cart

        def dire(self):

            dir = np.zeros((misil.t_impacto(self)[2]-1,3))

            dirx = np.zeros((misil.t_impacto(self)[2]-1,1))
            diry = np.zeros((misil.t_impacto(self)[2]-1,1))
            dirz = np.zeros((misil.t_impacto(self)[2]-1,1))

            n_frames = misil.t_impacto(self)

            coord_cart = misil.coord_misil(self)[1]
            # no estaba puesto el 1
            # tambien seria valido: a, coord_cart = misil.coord_misil(self)

            for i in range(misil.t_impacto(self)[2]-1):
                dir[i,0] = coord_cart[i+1,0]-coord_cart[i,0]
                dir[i,1] = coord_cart[i+1,1]-coord_cart[i,1]
                dir[i,2] = coord_cart[i+1,2]-coord_cart[i,2]
                
                if (i < n_frames[1]):
                    pass
                else:
                    dirz[i] = np.arccos(dir[i,2]/np.linalg.norm(dir[i]))
                    dirx[i] = np.arccos(dir[i,0]/(np.sin(dirz[i])*np.linalg.norm(dir[i])))
                    diry[i] = np.pi/2 - dirx[i]

            return dirx, diry, dirz 

    prueba1 = posInicial(x0,y0,z0,v0,100,1)
    pos_anim, r_meteo, pos_impacto = prueba1.vect_pos()
    prueba2 = misil(1,prueba1.tiempos_frames(), r_meteo, prueba1.coord() )

    message['pos_meteo'] = pos_anim.tolist()
    message['pos_misil']= prueba2.coord_misil()[0].tolist()
    message['dir_x']= prueba2.dire()[0].tolist()
    message['dir_y']= prueba2.dire()[1].tolist()
    message['dir_z']= prueba2.dire()[2].tolist()
    message['impacto'] = pos_impacto.tolist()

    print(message['pos_meteo'])
    print(message['pos_misil'])
    print(message['dir_x'])
    print(message['dir_y'])
    print(message['dir_z'])
    print(message['impacto'])

    json.dumps({'pos_meteo': message['pos_meteo']})
    json.dumps({'pos_misil': message['pos_misil']})
    json.dumps({'dir_x': message['dir_x']})
    json.dumps({'dir_y': message['dir_y']})
    json.dumps({'dir_z': message['dir_z']})
    json.dumps({'impacto': message['impacto']})

    emit('my_response',
         {'pos_meteo': message['pos_meteo'], 'pos_misil': message['pos_misil'], 'dir_x': message['dir_x'], 'dir_y': message['dir_y'], 'dir_z': message['dir_z'], 'impacto': message['impacto']})


if __name__ == '__main__':
    socketio.run(app)