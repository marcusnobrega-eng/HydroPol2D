function [output_matrix] = resize_matrix(input_matrix,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,dimension)
   % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    %                                                                 %
    %                 Produced by Marcus Nobrega Gomes Junior         %
    %                 e-mail:marcusnobrega.engcivil@gmail.com         %
    %                           September 2021                        %
    %                                                                 %
    %                 Last Updated: 11 September, 2021                %
    %                                                                 %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
%
%   -----  SYNTAX  -----
%   [output_matrix] = resize_matrix(input_matrix,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix);
%   -----  DESCRIPTION  -----
%   This function creates a new expanded output_matrix from an input_matrix
%   using a new coordinate system
%
%   -----  INPUT  -----
%  input_matrix(mandatory) : the MatLab matrix to be expanded in a new domain
%  
%  domain(mandatory) : expanded domain (zero matrix)
%  The other inputs are calculated using formulas in the script

%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
scratch_matrix = input_matrix(ymin_domain_abstracted:ymax_domain_abstracted,xmin_domain_abstracted:xmax_domain_abstracted,1:dimension);
output_matrix = repmat(domain,[1 1 dimension]);
output_matrix((dymin_matrix):(dymax_matrix),(dxmin_matrix):(dxmax_matrix),:) = scratch_matrix;
end

