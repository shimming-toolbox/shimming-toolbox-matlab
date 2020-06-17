%% HEADER
% This function crawls across the various folders of the project, and 
% create .md files. Each file contains the formatted documentation for 
% each function that is located under the folders.
% 
% Authors: Ryan Topfer, Alexandre D'Astous, Julien Cohen-Adad

%% Initial setup
% Add various paths
cd ..
addpath(genpath('./helpDocMd/src'))
addpath(genpath('./s'))

%% API doc
% overwrite shimming-toolbox/docs/contributing/api_documentation/
!mkdir s/docs/3_contributing/api_documentation

% Generate API documentation
src = './s';

%% Run documentation of API
src = Documentor.findfiles( src );

% Remove files not working
src = src(~contains(src,'s/Coils/Shim_Siemens/Shim_Prisma/ShimOpt_Prisma.m'));
src = src(~contains(src,'s/Coils/Shim_Siemens/Shim_Prisma/Shim_HGM_Prisma/ShimOpt_HGM_Prisma.m'));
src = src(~contains(src,'s/Coils/Shim_Siemens/Shim_Prisma/Shim_IUGM_Prisma_fit/ShimOpt_IUGM_Prisma_fit.m'));
src = src(~contains(src,'s/Ui/ShimUse.m'));

% Remove files not to be documented
src = src(~contains(src,'s/tests'));

% Call documentor
Options.outputDir = './s/docs/contributing/api_documentation';
Dr = Documentor( src , Options ) ;

% Generate documentation
Dr.printdoc() ;

disp('done')
exit;
                       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General purpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Read Yaml file
% filepath = './shimming-toolbox/mkdocs.yml';
% configFile = yaml.ReadYaml(filepath);
% %% Write Yaml file
% % yaml stores the config file a a data structure, to do so, it needs to
% % remove the spacings from certain nav fields. The following function
% % finds those fields and adds the spaces back
% yaml.WriteYaml('./shimming-toolbox/mkdocs2.yml', configFile);
% keywords = {'Getting Started'; 'API Documentation'; 'Other Ressources'};
% correctSpacings('./shimming-toolbox/mkdocs2.yml', keywords)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% function correctSpacings(filepath, keywords)
%     
% Takes a `filepath` to a ".yml" file and add spaces to the corresponding
% `keywords` when found in the ".yml" file
%
%     % create cell array containing keyword without spaces
%     for i = 1:numel(keywords)
%         noSpacing{i,1} = keywords{i}(~isspace(keywords{i}));
%     end
%     
%     % Read txt into cell A
%     fileID = fopen(filepath,'r');
%     i = 1;
%     while ~feof(fileID)
%         tline = fgetl(fileID);
%         
%         % In a line, look for each keywords
%         for iKeyword = 1:numel(keywords)
%              if(~isempty(strfind(tline,noSpacing{iKeyword})))
%                  
%                  keywordStarts = strfind(tline,noSpacing{iKeyword}) - 1;
%                  
%                  iSpace = 0;
%                  newKeyword = keywords{iKeyword};
%                  spaceInKeyword = strfind(newKeyword,' ') - 1 + iSpace;
%                  while(~isempty(spaceInKeyword))
%                      spaceInTline = keywordStarts + spaceInKeyword ;
%                      tline = [tline(1:spaceInTline), ' ' ,  tline(spaceInTline+1:end)] ;
%                      
%                      newKeyword = [newKeyword(1:spaceInKeyword-iSpace),  newKeyword(spaceInKeyword+2-iSpace:end)]
%                      iSpace = iSpace+1;
%                      spaceInKeyword = strfind(newKeyword,' ') - 1 + iSpace;
%                      
%                  end
%              end
%         end
% 
%         newText{i} = tline;
%         i = i+1;
%     end
%     fclose(fileID);
%     
%     % Write newText into txt
%     fileID = fopen(filepath, 'w');
%     for i = 1:numel(newText)
%         if i == numel(newText)
%             fprintf(fileID,'%s', newText{i});
%             break
%         else
%             fprintf(fileID,'%s\n', newText{i});
%         end
%     end
%     fclose(fileID);
% 
% end
