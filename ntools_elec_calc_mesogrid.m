function elec = ntools_elec_calc_mesogrid(elec,subjectpath)

% calculate grid 64-128 for PMT model 2110-128-021

mg = strncmpi('MG',elec(:,5),1);

if any(mg)
    elec_regular_grid = elec(~mg,:);
    elec = elec(mg,:);
    
    name = regexpi(elec(:,1),'[A-Za-z]*[^\d*]','match');
    name = unique([name{:}]','stable');
    
    elec_mg = cell(length(name),5); l = 0;
    
    for i=1:length(name)
        tf = strncmpi(name{i},elec(:,1),length(name{i}));
        % G01~G64
        elec_pos_mg = cell2mat(elec(tf,2:4));
        
        % determine the hemisphere that grid locates
        if elec_pos_mg(:,1)>0
            sph = 'rh';
        elseif elec_pos_mg(:,1)<0
            sph = 'lh';
        else
            error(['Grid looks across the hemisphere ', name{i}]);
        end
        
        %%
        % interpolate to 15*15
        [r,c] = meshgrid(1:8,1:8);
        [rr,cc] = meshgrid(1:0.5:8,1:0.5:8);
        X = interp2(r,c,reshape(elec_pos_mg(:,1),8,8)',rr,cc,'spline');
        Y = interp2(r,c,reshape(elec_pos_mg(:,2),8,8)',rr,cc,'spline');
        Z = interp2(r,c,reshape(elec_pos_mg(:,3),8,8)',rr,cc,'spline');
        
        % find G65~G128
        ind = false(15,15);
        ind(1,[2,4,6]) = 1; % g65~g67
        ind(2,1:6) = 1; % g68~g73
        ind(3,[2,4,6]) = 1; % g74~g76
        ind(4,1:9) = 1; % g77~g85
        ind(5,2:2:12) = 1; % g86~g91
        ind(6,1:14) = 1; % g92~g105
        ind(7,2:2:12) = 1; % g106~g111
        ind(9,8:2:14) = 1; % g112~g115
        ind(10,7:14) = 1; % g116~g123
        ind(11,8:2:14) = 1; % g124~g127
        ind(12,8) = 1; % g128

        xx = X'; xx = xx(:);
        yy = Y'; yy = yy(:);
        zz = Z'; zz = zz(:);
        ind = ind'; ind = ind(:);
        xx = xx(ind); yy = yy(ind); zz = zz(ind);
%         plot3(xx,yy,zz,'ro'); axis tight; axis equal;

        surf = fs_read_surf(fullfile(subjectpath,'surf', [sph '.pial']));
        if ~isfield(surf,'coords')
            surf.coords = surf.vertices;
        end
        kk = dsearchn(surf.coords,[xx,yy,zz]);
        elec_pos_mg = [elec_pos_mg;surf.coords(kk,:)];

        for j = 1:size(elec_pos_mg,1)
            elec_mg(l+j,1) = cellstr(sprintf('%s%.3d',char(name{i}),j));
            elec_mg(l+j,2:4) = num2cell(elec_pos_mg(j,:));
            if j<=64 
                elec_mg(l+j,5) = {'G'};
            else
                elec_mg(l+j,5) = {'EG'};
            end
        end
        
        l = j;

    end
    
    elec = [elec_mg;elec_regular_grid];

end

%% map EG to the surface

eg = strncmpi('EG',elec(:,5),1);
if any(eg)
    elec_others = elec(~eg,:);
    elec = elec(eg,:);
    
    name = regexpi(elec(:,1),'[A-Za-z]*[^\d*]','match');
    name = unique([name{:}]','stable');
    
    elec_eg = cell(length(name),5); l = 0;
    
    for i=1:length(name)
        tf = strncmpi(name{i},elec(:,1),length(name{i}));

        elec_pos_eg = cell2mat(elec(tf,2:4));
        
        % determine the hemisphere that grid locates
        if elec_pos_eg(:,1)>0
            sph = 'rh';
        elseif elec_pos_eg(:,1)<0
            sph = 'lh';
        else
            error(['Grid looks across the hemisphere ', name{i}]);
        end
        
        surf = fs_read_surf(fullfile(subjectpath,'surf', [sph '.pial']));
        if ~isfield(surf,'coords')
            surf.coords = surf.vertices;
        end
        kk = dsearchn(surf.coords,[elec_pos_eg(:,1),elec_pos_eg(:,2),elec_pos_eg(:,3)]);
        elec_pos_eg = surf.coords(kk,:);

        for j = 1:size(elec_pos_eg,1)
            if any(mg)
                elec_eg(l+j,1) = cellstr(sprintf('%s%.3d',char(name{i}),j+64)); % EG index: 065-128 
            else
                elec_eg(l+j,1) = cellstr(sprintf('%s%.3d',char(name{i}),j)); % EG index: 001-128
            end
            elec_eg(l+j,2:4) = num2cell(elec_pos_eg(j,:));
            elec_eg(l+j,5) = {'EG'};
        end
        
        l = j;

    end
    
    elec = [elec_others; elec_eg];

end
