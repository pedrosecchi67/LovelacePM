function main()
p1=[[0.02164787 0.36943459 0.79072794 0.88529385]
 [0.97427122 0.42364793 0.46298053 0.72752047]
 [0.18750825 0.28440592 0.09620873 0.59595417]]
p2=[[0.30597923 0.20344246 0.17598331 0.43976118]
 [0.10262271 0.96176601 0.95225361 0.32966035]
 [0.30095085 0.27343563 0.21355231 0.45733722]]
cp1=baricenter(p1);
cp2=baricenter(p2);
AIC1=RefQuadAIC(p1, cp1)
AIC2=RefQuadAIC(p2, cp1)
AIC3=RefQuadAIC(p1, cp2)
AIC4=RefQuadAIC(p2, cp2)
end
function b=baricenter(p)
b=[mean(p(1, :)); mean(p(2, :)); mean(p(3, :))];
end
function AIC=RefQuadAIC(p, ref)
pts=[p(:, 1)-ref, p(:, 2)-ref, p(:, 3)-ref, p(:, 4)-ref];
AIC=QuadAIC(pts);
end
function AIC=QuadAIC(p)
fprintf('Point coordinates:\n');
p
AIC=LineIntegral(p(:, 1), p(:, 2))+LineIntegral(p(:, 2), p(:, 3))+LineIntegral(p(:, 3), p(:, 4))+LineIntegral(p(:, 4), p(:, 1));
end
function vbar=LineIntegral(p1, p2)
fprintf('=================\n');
fprintf('Line Integral between [%f %f %f] and [%f %f %f]\n', [p1; p2]);
vbar=zeros(3, 1);
u=p2-p1;
th1=ang(p1, u);
th2=ang(p2, u);
R=norm(p1)*cos(th1);
vers=cross(p1, u);
vers=vers/norm(vers);
fprintf('thetas: %f %f\n', [th1, th2]);
fprintf('R: %f\n', [R]);
vbar=-vers*sin(th1)/(4*pi*R);
vers=cross(p2, u);
vers=vers/norm(vers);
vbar=vbar+vers*sin(th2)/(4*R*pi);
fprintf('Total disturbance velocity: %f %f %f\n', vbar);
fprintf('=================\n');
end
function alpha=ang(u, v)
alpha=asin(u'*v/(norm(u)*norm(v)));
end
