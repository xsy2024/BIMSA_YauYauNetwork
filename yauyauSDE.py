from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize

        
class yauyauSDE:
    
    
    
    
    
    def __init__(self, name, data):
        """
        Initialize the yauyauSDE object with dataset and names.
        
        Parameters:
            name: list or other iterable containing names for the datasets.
            data: dictionary or similar structure containing the data.
        """
       
        self.Total_name = name           # List of dataset names
        self.Total_data = data           # Dictionary containing the data for each name
        self.Total_weight = dict()       # Dictionary to store calculated weights
        self.Total_weight_ = dict()      # Additional weight storage
        self.static_w = dict()           # Dictionary for static weight (accumulated weights)
        self.intercept_self_ = dict()    # Dictionary for storing self intercepts
        self.coef_ = dict()              # Dictionary for regression coefficients
        self.intercept_ = dict()         # Dictionary for regression intercepts
        self.R2 = dict()                 # Dictionary to store R^2 scores for regressions
        self.gene = dict()               # Dictionary to store predicted gene values
        self.gene_self = dict()          # Dictionary for self gene predictions
        self.mse = dict()                # Dictionary for storing mean squared errors
        self.W = dict()                  # Dictionary for weight matrices computed by different methods
        self.W_p = dict()                # Dictionary for predicted weights
        
    
        
        
    def LegendrePolynomials(self,N,x):
       """
        Compute the Legendre polynomial of order N evaluated at x.
        
        Parameters:
            N: int, the order of the Legendre polynomial.
            x: numeric or symbolic variable at which the polynomial is evaluated.
        
        Returns:
            The value of the Legendre polynomial of degree N at x.
        """
        if N == 0:
            return 1
        if N == 1:
            return x
        p0 = self.LegendrePolynomials(0,x)
        p1 = self.LegendrePolynomials(1,x)
        assert N>=2
        for i in range(1,N):
            p = (2*i+1)/(i+1)*x*p1 - i/(i+1)*p0
            p0 = p1
            p1 = p
        return p1   
    
    def sub(self):
        """
        Preprocess the Total_data for each dataset name by subtracting the value at time 0
        and selecting the first 17 columns.
        """
        for j in self.Total_name:
            # Subtract the first row (at time 0) from each row along axis 0
            self.Total_data[j] = self.Total_data[j].sub(self.Total_data[j][0.0], axis=0)
            self.Total_data[j] = self.Total_data[j].iloc[:,:17]
  
    def L_conti(self, no,step): 
        """
        Fit a polynomial (up to 4th degree) to the data for a given sample and generate
        new predictions with a given step size.
        
        Parameters:
            no: index of the sample in the data.
            step: float, step size for new time grid.
        """
        t = np.arange(0., 4.0001, 0.25)
        t_2 = t**2
        t_3 = t**3
        t_4 = t**4
        X = np.vstack((t, t_2, t_3, t_4))
        X = X.T
        t_new = np.arange(0., 4.0001, step)
        X_new = np.vstack((t_new, t_new**2, t_new**3, t_new**4))
        X_new = X_new.T
        gene = []
        for i in range(3):
            y = np.array(self.Total_data[self.Total_name[i]].iloc[no,:])
            model = LinearRegression(fit_intercept=False)
            model.fit(X,y)           
            y_pred = model.predict(X_new)
            gene.append(y_pred)
        self.gene[no] = gene
        
    def h(self):
        """
        Compute mean values and additional features from the Total_data,
        fit a regression model and calculate weights.
        """
        X_0 = np.array([self.Total_data[j].mean(0) for j in self.Total_name])
       
        X_1 = np.array([self.Total_data[j].mean(0)**2 for j in self.Total_name])
        X_2 = np.array([self.Total_data[j].mean(0)**3 for j in self.Total_name])
        X_3 = np.array([np.sin(self.Total_data[j].mean(0)) for j in self.Total_name])
        X_4 = np.array([np.cos(self.Total_data[j].mean(0)) for j in self.Total_name])
        X_5 = np.array([np.cos(2*self.Total_data[j].mean(0)) for j in self.Total_name])
        X_6 = np.array([np.cos(3*self.Total_data[j].mean(0)) for j in self.Total_name])
        X_7 = np.array([np.sin(2*self.Total_data[j].mean(0)) for j in self.Total_name])
        X = np.concatenate((X_0, X_1, X_2), axis = 0)
        X = np.array(X).T
        print(X.shape)
        y = np.array([self.Total_data[f"d_{j}"].mean(0).iloc[:17] for j in self.Total_name])
        
        coef = []
        intercept = []
        y_preds = []
        for i in range(3):
             model = LinearRegression()
             model.fit(X, y[i])
             y_pred = model.predict(X)
             mse = mean_squared_error(y[i], y_pred)
             r2 = r2_score(y[i], y_pred)
             residuals = y[i] - y_pred
             #print(mse)
             plt.plot(y[i])
             print(mse)
             print(r2)
             coef.append(model.coef_)
             intercept.append(model.intercept_)
             y_preds.append(y_pred)
             
        self.coef = coef
        self.intercept = intercept
        self.y_preds = y_preds
        
        w = np.zeros((17,3,3))
        for i in range(3):
            for j in range(3):               
                w[:,i,j] = sum([X[:,3*j+m]*coef[i][3*j+m] for m in range(3)])
                                   
        self.Total_weight["mean"] = w
        
        
        np.savetxt("coef.csv",self.coef,delimiter=",") 
        np.savetxt("intercept.csv",self.coef,delimiter=",") 
        np.savetxt("Total_data.csv",X_0,delimiter=",")
    
    def w_individual(self,state,name,no):
        """
        Compute individual weight matrices based on the provided state.
        
        Parameters:
            state: np.array, the state variable vector.
            name: str, the key under which to store computed weights.
            no: sample index.
        """
        X = np.concatenate((state,state**2,np.cos(state)),axis = 1)
        #X = np.concatenate((self.LegendrePolynomials(1,state),self.LegendrePolynomials(2,state),self.LegendrePolynomials(3,state)),axis = 1)
        length = state.shape[0]
        w_s = np.zeros((length,3,3))
        for i in range(3):
            for j in range(3):               
                w_s[:,i,j] = sum([X[:,3*m+j]*self.coef_[no][i][3*m+j] for m in range(3)])
        self.Total_weight[name] = w_s
        esp = ["E","P","S"]
        color = [(0.1,0.1,0.8,0.8),(0.8,0.1,0.1,0.8), (0.1,0.8,0.1,0.8)]
        for i in range(3):
            plt.figure() 
            for j in range(3):
                plt.plot(np.arange(0,4,4/length),w_s[:,i,j], label = f"{esp[j]} to {esp[i]} {no}", color = color[j] )
            np.savetxt(f'{name}_{i}.csv', w_s[:,i,:], delimiter = ",")
            plt.xlabel("Time")
            plt.ylabel("Interaction Strength")
            plt.legend()
            plt.savefig(f"{no}_{esp[i]}")
            
        step = 4/length
        self.static_w[name] = np.sum(w_s,axis=2)*step
        
    def w_deter(self,name,no):
        """
        Determine weight matrices using a deterministic approach based on gene predictions.
        
        Parameters:
            name: str, key for storing weight matrices.
            no: sample index.
        """
        #state = np.array([self.Total_data[j].iloc[no,:17] for j in self.Total_name]).T
        state = np.vstack(self.gene[no]).T
        print(state.shape)
        #X = np.concatenate((state,state**2,np.cos(state)),axis = 1)
        X = np.concatenate((state,state**2,state**3),axis = 1)
        length = state.shape[0]
        w_s = np.zeros((length,3,3))
        for i in range(3):
            for j in range(3):               
                w_s[:,i,j] = sum([X[:,3*m+j]*self.coef_[no][i][3*m+j] for m in range(3)])
        self.Total_weight[name] = w_s
        esp = ["E","P","S"]
        color = [(0.1,0.1,0.8,0.8),(0.8,0.1,0.1,0.8), (0.1,0.8,0.1,0.8)]
        for i in range(3):
            plt.figure() 
            for j in range(3):
                plt.plot(np.arange(0,4,4/length),w_s[:,i,j], label = f"{esp[j]} to {esp[i]} {no}", color = color[j] )
            
            np.savetxt(f'{name}_{i}.csv', w_s[:,i,:], delimiter = ",")
            plt.xlabel("Time")
            plt.ylabel("Interaction Strength")
            plt.legend()
            plt.savefig(f"deter_{no}_{esp[i]}")
        
        step = 4/length
        self.static_w[name] = np.sum(w_s,axis=0)*step
        
    def w(self,state,name):
        """
        Compute weight matrices using mean regression coefficients.
        
        Parameters:
            state: np.array, state variable.
            name: str, key for storing computed weights.
        """
        X = np.concatenate((state,state**2,np.cos(state)),axis = 1)
        length = state.shape[0]
        w_s = np.zeros((length,3,3))
        for i in range(3):
            for j in range(3):               
                w_s[:,i,j] = sum([X[:,3*j+m]*self.coef[i][3*j+m] for m in range(3)])
        self.Total_weight[name] = w_s
        esp = ["E","P","S"]
        color = [(0.1,0.1,0.8,0.8),(0.8,0.1,0.1,0.8), (0.1,0.8,0.1,0.8)]
        for i in range(3):
            plt.figure() 
            for j in range(3):
                plt.plot(np.arange(0,4,4/length),w_s[:,i,j], label = f"{esp[j]} to {esp[i]} mean", color = color[j] )
            
            plt.legend()
            plt.savefig(f"mean_{esp[i]}")
            
        step = 4/length
        self.static_w_ = np.sum(w_s,axis=0)*step
    
    def w_new(self, state, no, name):
        """
        Compute new weight matrices using gene predictions and adjust for noise.
        
        Parameters:
            state: np.array, the original state variable.
            no: sample index.
            name: str, key for storing computed weights.
        """
        state_ = np.vstack(self.gene[no]).T
        X = np.concatenate((state,state**2,state**3),axis = 1)
        X_ = np.concatenate((state_,state_**2,state_**3),axis = 1)
        length = X.shape[0]
        w_s = np.zeros((length,3,3))
        w_s_ = np.zeros((length,3,3))
        step = 4/length
        A = [[i+j*3 for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):   
                if j == 2:
                    w_s[:,i,j] = 0.45*np.dot(X[:,A[j]],self.coef_[no][i][A[j]])
                else:
                    w_s[:,i,j] = np.dot(X[:,A[j]],self.coef_[no][i][A[j]])
                w_s_[:,i,j] = np.dot(X_[:,A[j]],self.coef_[no][i][A[j]])
            
            np.savetxt(f"{name}_{i}.csv", w_s_[:,i,:], delimiter = ",")    
            np.savetxt(f"{name}_{i}_noise.csv", w_s[:,i,:], delimiter = ",")
            
        self.Total_weight[name] = w_s
        self.Total_weight_[name] = w_s_
        
        esp = ["E","P","S"]
        color = [(0.1,0.1,0.8,0.8),(0.8,0.1,0.1,0.8), (0.1,0.8,0.1,0.8)]
        for i in range(3):
            plt.figure() 
            for j in range(3):
                plt.plot(np.arange(0,4,4/length),w_s[:,i,j], label = f"{esp[j]} to {esp[i]}", color = color[j] )
                plt.plot(np.arange(0,4,4/length),w_s_[:,i,j], label = f"{esp[j]} to {esp[i]}", color = color[j] )

            plt.legend()
            plt.savefig(f"new_{esp[i]}")
            
        #self.static_w_ = np.sum(w_s,axis=0)*step
    
    def H(self, no):
        """
        Fit a regression model for the sample with index 'no' using various feature transformations.
        
        Parameters:
            no: sample index.
        """
        X_0 = np.array([self.Total_data[j].iloc[no,:] for j in self.Total_name])
        X_1 = np.array([self.Total_data[j].iloc[no,:]**2 for j in self.Total_name])
        X_2 = np.array([self.Total_data[j].iloc[no,:]**3 for j in self.Total_name])
        X_3 = np.array([np.sin(self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X_4 = np.array([np.cos(self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        
        X_5 = np.array([self.Total_data[j].iloc[no,:]**4 for j in self.Total_name])
        X_6 = np.array([np.sin(2*self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X_7 = np.array([np.cos(2*self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X = np.concatenate((X_0, X_1, X_2), axis = 0)
        state = np.vstack(self.gene[no])
        X = np.concatenate((state,state**2,state**3),axis = 0).T
        y = np.array([self.Total_data[f"d_{j}"].iloc[no,:17] for j in self.Total_name for j in self.Total_name])
        
        coef = []
        intercept = []
        R2 = []
        for i in range(3):
             model = LinearRegression()
             model.fit(X, y[i])
             y_pred = model.predict(X)
             mse = mean_squared_error(y[i], y_pred)
             r2 = r2_score(y[i], y_pred)
             residuals = y[i] - y_pred
             #print(mse)
             plt.plot(y[i])
             #print(r2)
             #print(mse)
             coef.append(model.coef_)
             intercept.append(model.intercept_)
             R2.append(r2)
             
        self.coef_[no] = coef
        self.intercept_[no] = intercept
        self.R2[no] = R2
        
    def H_self(self, no):
        """
        Similar to H(), but uses 'self' variant for internal gene predictions.
        
        Parameters:
            no: sample index.
        """
        X_0 = np.array([self.Total_data[j].iloc[no,:] for j in self.Total_name])
        X_1 = np.array([self.Total_data[j].iloc[no,:]**2 for j in self.Total_name])
        X_2 = np.array([self.Total_data[j].iloc[no,:]**3 for j in self.Total_name])
        X_3 = np.array([np.sin(self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X_4 = np.array([np.cos(self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        
        X_5 = np.array([self.Total_data[j].iloc[no,:]**4 for j in self.Total_name])
        X_6 = np.array([np.sin(2*self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X_7 = np.array([np.cos(2*self.Total_data[j].iloc[no,:]) for j in self.Total_name])
        X = np.concatenate((X_0, X_1, X_2), axis = 0)
        X = np.array(X).T
        state = np.vstack(self.gene[no])
        X = np.concatenate((state,state**2,state**3),axis = 0).T
        length = X.shape[0]
        y = np.array([self.Total_data[f"d_{j}"].iloc[no,:17] for j in self.Total_name for j in self.Total_name])
        
        coef = []
        intercept1 = []
        intercept2 = []
        R2 = []
        residuals = []
        gene = []
        A = [[i+j*3 for j in range(3)] for i in range(3)]
        W = []
        W_p = []
        Total = []
        name = ["E","P","S"]
        cr = ["purple","orange",'g']
        for i in range(3):
             c = np.zeros(9)
             w_s= np.zeros((length,3))
             w_p= np.zeros((length,3))
             number = np.array(range(9))
             
             a = [i+j*3 for j in range(2)]
             ai = np.delete(number, A[i])
             X_self = X[:,a]
             X_dep = np.delete(X,A[i],axis = 1)
                       
             model = LinearRegression(fit_intercept=False)
             model.fit(X_self, y[i])
             
             c[a] = model.coef_
             y_pred = model.predict(X_self)
             
             mse = mean_squared_error(y[i], y_pred)
             t_eval = np.linspace(0,4,1700)
             t_span = (0,4)
             y0 = [0]
             y_hat = [2*(4/length)*sum(y_pred[:t]) for t in range(length)]
             model = LinearRegression(fit_intercept=False)
             model.fit(X_self,y_hat)
             c[a] = model.coef_
             gene.append(y_hat)
             plt.figure()
             plt.title(name[i])
             plt.plot(np.linspace(0,4,length),y_hat[:length], label = "Indep_effect", color = "r")
             plt.plot(np.linspace(0,4,length),np.zeros(length),linestyle = "dashed",color = "k")
             plt.scatter(np.linspace(0,4,length), X_0[i,:length], color = "k", s = 20, alpha=0.8)
             
             w_s[:,i] = y_hat[:length]
             w_p[:,i] = y_hat[:length]
             y_obs = X[:,i] - y_hat
             y_obs = y_obs[:length]
             X_dep = X_dep[:length,:]
             intercept1.append(model.intercept_)
             model = LinearRegression(fit_intercept=False)
             model.fit(X_dep, y_obs)
             y_pred = model.predict(X_dep)
             plt.plot(np.linspace(0,4,length),w_s[:,i]+y_pred, color = "b")
             Total.append(w_s[:,i]+y_pred)
             
             c[ai] = model.coef_
             print(f"model.coef_ : {model.coef_}")
             print(f"c:{c}")
             
             for m in range(3):
                 cm = np.zeros(9)
                 model = LinearRegression(fit_intercept=False)
                 if m != i:
                     cm[A[m]] = c[A[m]]
                     print(f"cm:{cm}")
                     print(f"c:{c}")
                     X_dep = X[:17,A[m]]
                     
                     #print(X_dep.shape)
                     model.fit(X_dep, y_obs)
                     y_pred = model.predict(X_dep)
                     w_p[:,m] = y_pred
                     c[A[m]] = model.coef_
                     #plt.plot(y_pred)
                     y_obs = y_obs-y_pred
                     dm = np.dot(X[:length,:], cm)
                     #w_s[:,m] = [0.25*sum(dm[:t]) for t in range(17)]
                     w_s[:,m] = dm
                     plt.plot(np.linspace(0,4,length),w_p[:,m], label = name[m], color = cr[m])
                     #plt.plot(w_s[:,m], label = name[m])
             plt.legend()
             
             mse = mean_squared_error(X[:length,i], y_pred+y_hat[:length])
             np.savetxt(f"{no}_W_{i}.csv",w_p,delimiter = ",")
             W.append(w_s)
             W_p.append(w_p)
             
             coef.append(c)
             
        
        self.gene_self[no] = np.array(gene)    
        self.coef_[no] = coef
        self.mse[no] = mse
        self.W[no] = W
        self.W_p[no] = W_p
        
        np.savetxt(f"{no}_t_fit.csv",state,delimiter = ",") 
        np.savetxt(f"{no}_c.csv",coef,delimiter = ",")
        np.savetxt(f"{no}_total.csv",Total,delimiter = ",")
        
    def ode(self, y, t, p):
        """
        ODE function used in integration.
        
        Parameters:
            y: current value of dependent variable.
            t: time variable.
            p: parameter vector [a1, a2, a3].
        
        Returns:
            dydt: derivative computed as a1*y + a2*y^2 + a3*y^3.
        """
        a1 ,a2, a3 = p
        dydt = a1*y+a2*y**2 +a3*y**3
        
        #dydt = a1*y*(1-y/(a2+1e-5))**a3
        print(dydt)
        return dydt
    
    def odeerror(self, p, y0, t_span, t_eval, y_obs):
        """
        Compute error between ODE solution and observed data.
        
        Parameters:
            p: parameter vector for ODE.
            y0: initial condition.
            t_span: tuple (t0, tf).
            t_eval: time points where solution is evaluated.
            y_obs: observed data.
        
        Returns:
            error: sum of squared differences between solution and observed data.
        """
        solve = solve_ivp(self.ode,  t_span, y0, "RK45",t_eval, args=(p,))
        y_hat = solve.y[0][0:-1:100]
        error = np.sum((y_hat-y_obs)**2)
        return error
        
        
    def f(self, x, no):
        """
        Compute the polynomial function based on stored coefficients and intercept for sample no.
        
        Parameters:
            x: input feature vector.
            no: sample index.
        
        Returns:
            The computed value of the polynomial function.
        """
        print(x**2)
        return np.array(self.intercept_[no])+ np.dot(np.array(self.coef_[no])[:,:3],x) + np.dot(np.array(self.coef_[no])[:,3:6],x**2) +np.dot(np.array(self.coef_[no])[:,6:9],x**3)
    
    def f_mean(self, x):
        """
        Compute the mean polynomial function using mean coefficients and intercept.
        
        Parameters:
            x: input feature vector.
        
        Returns:
            Computed function value.
        """
        return np.array(self.intercept)+ np.dot(np.array(self.coef)[:,:3],x) + np.dot(np.array(self.coef)[:,3:6],x**2) + np.dot(np.array(self.coef)[:,6:],np.cos(x))
    
    def Euler_gene(self, t, dt, no, x_0 = 0):
        """
        Solve an ODE using Euler's method to predict gene expression.
        
        Parameters:
            t: total time.
            dt: time step.
            no: sample index.
            x_0: initial condition (default=0).
        """
        step = int(t/dt)
        x = np.zeros((3,step+1))
        for j in range(step):
            x[:,j+1] = x[:,j] +self.f(x[:,j], no)*dt
            
        self.gene[no] = x
        name = ["E","P","S"]
        for i in range(3):
            plt.figure()
            plt.plot(np.arange(0,4.25,0.25),self.Total_data[self.Total_name[i]].iloc[no,:17], label = f"{no} real data {name[i]}", color = "black")
            plt.plot(np.arange(0,4.0001,0.001),x[i,:], label = "Fitting result", color = "blue")
            
            plt.xlabel("Time")
            plt.ylabel("Abundance")
            plt.legend()
            plt.savefig(f"Fitting_{no}_{name[i]}")
        np.savetxt(f"{no}_euler.csv", x,delimiter = ",")    
        np.savetxt(f'{no}_data.csv',np.array([self.Total_data[self.Total_name[i]].iloc[no,:17] for i in range(3)]),delimiter = "," )
    
    def Euler_gene_mean(self, t, dt, x_0 = 0):
        """
        Compute mean gene predictions using Euler's method with mean regression function.
        
        Parameters:
            t: total time.
            dt: time step.
            x_0: initial condition (default=0).
        """
        step = int(t/dt)
        x = np.zeros((3,step+1))
        for j in range(step):
            x[:,j+1] = x[:,j] +self.f_mean(x[:,j])*dt
            
        self.gene_ = x
        name = ["E","P","S"]
        for i in range(3):
            plt.figure()
            plt.plot(np.arange(0,4.25,0.25),self.Total_data[self.Total_name[i]].mean(0)[:4], label = f"real data {name[i]}", color = "black")
            plt.plot(np.arange(0,4.0001,0.001),x[i,:], label = "Sindy result", color = "blue")
            plt.legend()
            plt.savefig(f"Sindy_mean_{name[i]}")
    
