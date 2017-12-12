procedure periodsplinecoeffns_przed (n         : Integer;
                             interval_x,interval_f       : vectorprz;
                               var interval_a     : matrixprz;
                               var st    : Integer);
{---------------------------------------------------------------------------}
{                                                                           }
{  The procedure periodsplinecoeffns calculates the coefficients of the     }
{  periodic cubic spline interpolant for a function given by its values at  }
{  nodes.                                                                   }
{  Data:                                                                    }
{    n  - number of interpolation nodes minus 1 (the nodes are numbered     }
{         from 0 to n),                                                     }
{    x  - an array containing the values of nodes,                          }
{    f  - an array containing the values of function (the value of f[0]     }
{         should be equal to the value of f[n]).                            }
{  Result:                                                                  }
{    a  - an array of spline coefficients (the element a[k,i] contains the  }
{         coefficient before x^k, where k=0,1,2,3, for the interval         }
{         <x[i], x[i+1]>; i=0,1,...,n-1).                                   }
{  Other parameters:                                                        }
{    st - a variable which within the procedure periodsplinecoeffns is      }
{         assigned the value of:                                            }
{           1, if n<1,                                                      }
{           2, if there exist x[i] and x[j] (i<>j; i,j=0,1,...,n) such      }
{              that x[i]=x[j],                                              }
{           3, if f[0]<>f[n],                                               }
{           0, otherwise.                                                   }
{         Note: If st<>0, then the elements of array a are not calculated.  }
{  Unlocal identifiers:                                                     }
{    vector  - a type identifier of extended array [q0..qn], where q0<=0    }
{              and qn>=n,                                                   }
{    vector1 - a type identifier of extended array [q1..qn], where q1<=1    }
{              and qn>=n,                                                   }
{    vector2 - a type identifier of extended array [q1..qn1], where q1<=1   }
{              and qn1>=n-1,                                                }
{    vector3 - a type identifier of extended array [q2..qn], where q2<=2    }
{              and qn>=n,                                                   }
{    matrix  - a type identifier of extended array [0..3, q0..qn1], where   }
{              q0<=0 and qn1>=n-1.                                          }
{                                                                           }
{---------------------------------------------------------------------------}
var i,k        : Integer;
    w,v,y,z,xi : Interval;
    u          : array of Interval;
    b,c,d      : array of Interval;
    p          : array of Interval;
    q          : array of Interval;
begin

   SetLength(u,n+1);
   SetLength(b,n+1);
   SetLength(c,n+1);
   SetLength(d,n+1);
   SetLength(p,n+1);
   SetLength(q,n+1);
   SetLength(interval_a,4,n);

  if n<1
    then st:=1
    else if ((interval_f[0].a<>interval_f[n].a) OR (interval_f[0].b<>interval_f[n].b))
           then st:=3
           else begin
                  st:=0;
                  i:=-1;
                  repeat
                    i:=i+1;
                    for k:=i+1 to n do
                      if ((interval_x[i].a=interval_x[k].a) AND (interval_x[i].b=interval_x[k].b))
                        then st:=2
                  until (i=n-1) or (st=2)
                end;
  if st=0
    then begin
           if n>1
             then begin
                    v:=isub(interval_x[1],interval_x[0]);
                    z:=isub(interval_x[n],interval_x[n-1]);
                    b[n]:=idiv(v,iadd(z,v));
                    c[n]:=isub(int_read('1'),b[n]);
                    y:=interval_f[n];
                    d[n]:=idiv(imul(int_read('6'),isub(idiv(isub(interval_f[1],y),v),idiv(isub(y,interval_f[n-1]),z))),iadd(z,v));
                    for i:=1 to n-1 do
                      begin
                        z:=interval_x[i];
                        y:=isub(interval_x[i+1],z);
                        z:=isub(z,interval_x[i-1]);
                        v:=interval_f[i];
                        b[i]:=idiv(y,iadd(y,z));
                        c[i]:=isub(int_read('1'),b[i]);
                        d[i]:=idiv(imul(int_read('6'),isub(idiv(isub(interval_f[i+1],v),y),idiv(isub(v,interval_f[i-1]),z))),iadd(y,z))
                      end;
                    if n>2
                      then begin
                             u[1]:=int_read('2');
                             c[2]:=idiv(c[2],int_read('2'));
                             q[2]:=isub(int_read('0'),idiv(b[n],int_read('2')));
                             for i:=2 to n-2 do
                               begin
                                 v:=isub(int_read('2'),imul(b[i-1],c[i]));
                                 c[i+1]:=idiv(c[i+1],v);
                                 q[i+1]:=isub(int_read('0'),idiv(imul(q[i],b[i-1]),v));
                                 u[i]:=v
                               end;
                             v:=isub(int_read('2'),imul(c[n-1],b[n-2]));
                             q[n]:=idiv(isub(c[n],imul(q[n-1],b[n-2])),v);
                             u[n-1]:=v;
                             p[1]:=c[1];
                             for i:=2 to n-2 do
                               p[i]:=isub(int_read('0'),imul(c[i],p[i-1]));
                             p[n-1]:=isub(b[n-1],imul(c[n-1],p[n-2]));
                             v:=isub(int_read('2'),imul(c[1],p[2]));
                             for i:=2 to n-2 do
                               v:=isub(v,imul(p[i],p[i+1]));
                             u[n]:=isub(v,imul(p[n-1],q[n]));
                             for i:=2 to n-1 do
                               d[i]:=isub(d[i],imul(c[i],d[i-1]));
                             v:=d[n];
                             for i:=2 to n do
                               v:=isub(v,imul(q[i],d[i-1]));
                             d[n]:=v;
                             u[n]:=idiv(d[n],u[n]);
                             u[n-1]:=idiv(isub(d[n-1],imul(p[n-1],u[n])),u[n-1]);
                             for i:=n-2 downto 1 do
                               u[i]:=idiv(isub(isub(d[i],imul(b[i],u[i+1])),imul(p[i],u[n])),u[i]);
                           end
                      else begin
                             y:=d[1];
                             z:=d[2];
                             w:=isub(int_read('4'),imul(c[2],b[1]));
                             u[1]:=idiv(isub(imul(int_read('2'),y),imul(b[1],z)),w);
                             u[2]:=idiv(isub(imul(int_read('2'),z),imul(c[2],y)),w);
                           end
                  end
             else u[1]:=int_read('0');
           u[0]:=u[n];
           for i:=0 to n-1 do
             begin
               w:=interval_f[i];
               xi:=interval_x[i];
               z:=isub(interval_x[i+1],xi);
               y:=u[i];
               v:=isub(idiv(isub(interval_f[i+1],w),z),idiv(imul(iadd(imul(int_read('2'),y),u[i+1]),z),int_read('6')));
               z:=idiv(isub(u[i+1],y),imul(int_read('6'),z));
               y:=idiv(y,int_read('2'));
               interval_a[0,i]:=iadd(imul(isub(imul(iadd(imul(isub(int_read('0'),z),xi),y),xi),v),xi),w);
               w:=imul(imul(int_read('3'),z),xi);
               interval_a[1,i]:=iadd(imul(isub(w,imul(int_read('2'),y)),xi),v);
               interval_a[2,i]:=isub(y,w);
               interval_a[3,i]:=z
             end
         end
end;