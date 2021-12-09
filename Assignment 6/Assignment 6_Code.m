clc
k = 180;
model = createpde('thermal','steadystate');
importGeometry(model,'Triangular_fin.stl');

%for geometry model plot
figure(1);
pdegplot(model,'Facelabels','on','Facealpha',0.5);

%for Mesh model plot
msh = generateMesh(model,'Hmax',2);
figure(2);
pdeplot3D(model);

%specifying thermal properties and boundary conditions
thermalProperties(model, 'ThermalConductivity', k);
thermalBC(model,'Face',3,'Temperature',100);
thermalBC(model,'Face',[1 2 4 5],'ConvectionCoefficient',15,'AmbientTemperature',25);
Result = solve(model);

%for colormap 3D plot
figure(3);
pdeplot3D(model,'ColorMapData', Result.Temperature);