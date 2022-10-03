%Este programa halla la solución a la ecuación f(x)=0 usando el método de
%raíces múltiples

function r=RaicesMultp1(f,Xo,Nmax,tol)

gx=sym(f); %entrada simbólica de funcion
%Ingresar de la siguiente manera: f=@(Xsig) 'Función continua'
gx1=diff(gx); %primera derivada
gx2=diff(gx,2); %segunda derivada
i=1; %iteraciones

while i<=Nmax 
    Xsig=Xo-(f(Xo)*subs(gx1,Xo))/(((subs(gx1,Xo))^2)-f(Xo)*subs(gx2,Xo));%formula auxiliar Newton mejorada
    Error=abs(Xsig-Xo); %Error
    fprintf('%4.0f  %4.5f  %4.5f  %4.5f\n',i,Xsig,f(Xsig),Error); 
    if Error<tol
        disp('terminó, a continación se muestra el número de iteraciones y el valor Error');
        break;
    end
i=i+1;
Xo=Xsig;
    
end
r=fprintf('%4.0f,  %4.4f',i,Error); %Imprime al final el número de iteraciones y el valor de error en la iteración que se detuvo
disp('   Y el resultado de X:')
r=double(Xsig); %Valor de X donde se detuvo el programa
