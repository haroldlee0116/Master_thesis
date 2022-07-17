%% Figure4_7: similarity of the model between each validation
% download the 'recording3' data set from:
% https://drive.google.com/file/d/1IoR_z3wgqlxbCxOopUt_hs7aoR8j3_3I/view?usp=sharing
load recording3
V=recording3.V;
for i=1:7
    for j=1:7
        T_i = cpdgen([repmat(V(2*(i-1)+1),1,d),V{2*i}]);
        T_j = cpdgen([repmat(V(2*(j-1)+1),1,d),V{2*j}]);
        sim(i,j)= frob(T_i-T_j)/frob(T_i);
    end
end
h=heatmap(sim);
xlabel('validation realization');
ylabel('validation realization');



