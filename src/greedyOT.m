clc
clear all

load("rt40_mat_files\DOA_rt40.mat")
load("rt40_mat_files\pressure_rt40.mat")

length_signal = 19999;
start_signal = 206;
speed_of_sound = 343/48000;

pressure_cell = cell(23,2);
DOA_cell = cell(23,2);

pressure_cell{1,1} = P{1,1};
pressure_cell{1,2} = P{1,2};
DOA_cell{1,1} = DOA{1,1};
DOA_cell{1,2} = DOA{1,2};

for k=1:2
    for j=1:11
        DOA1 = DOA{j,k}(start_signal:start_signal+length_signal,:);
        P1 = P{j,k}(start_signal:start_signal+length_signal,:);
        DOA3 = DOA{j+1,k}(start_signal:start_signal+length_signal,:);
        P3 = P{j+1,k}(start_signal:start_signal+length_signal,:);
        
        for i=1:length(DOA1)
            DOA1_length(i,:) = DOA1(i,:)*(speed_of_sound*i);
            P1_adjusted(i,:) = P1(i,:)/(speed_of_sound*i);
        end
        for i=1:length(DOA3)
            DOA3_length(i,:) = DOA3(i,:)*(speed_of_sound*i);
            P3_adjusted(i,:) = P3(i,:)/(speed_of_sound*i);
        end
        
        DataCloud1 = [P1_adjusted, DOA1_length,P1,DOA1];
        DataCloud1 = sortrows(DataCloud1, -1);
        DataCloud3 = [P3_adjusted, DOA3_length,P3,DOA3];
        DataCloud3 = sortrows(DataCloud3, -1);
        
        pressureGreedy = (DataCloud1(:,5)+DataCloud3(:,5))/2;
        for i=1:length(DOA1)
            distanceGreedy(i,:) = vecnorm((DataCloud1(i,2:4)+DataCloud3(i,2:4))/2);
        end
        pointGreedy = ((DataCloud1(:,2:4)+DataCloud3(:,2:4))/2)./distanceGreedy;
        DataCloudGreedy = [pressureGreedy, distanceGreedy, pointGreedy];
        DataCloudGreedy = sortrows(DataCloudGreedy,2);
        pressure = DataCloudGreedy(:,1);
        DOA_greedy = DataCloudGreedy(:,3:5);
    
        pressure_cell{2*j,k} = pressure
        DOA_cell{2*j,k} = DOA_greedy
        pressure_cell{2*j+1,k} = P{j,1}
        DOA_cell{2*j+1,k} = DOA{j,1}
    
        clear DataCloud1
        clear DataCloud3
        clear pressureGreedy
        clear pointGreedy
        clear DataCloudGreedy
        clear pressure
        disp(j)
    end
end

save('interpolated_mat_files\P_interpolated.mat', "pressure_cell");
save("interpolated_mat_files\DOA_interpolated.mat", "DOA_cell");