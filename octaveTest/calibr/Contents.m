% Calibration toolbox v3.0  10-17-00
%
% Calibration:
%
% cacal    calibration routine for solving camera parameters
% cacalw   calibration routine using weighted least squares
% invmodel computes the parameters of the inverse distortion model
% imcorr   corrects distorted image coordinates
% imdist   adds radial and tangential distortion to image coordinates
% extcal   computes the extrinsic camera parameters for a calibrated camera
% configc  gives the configuration information
%
% Subroutines:
%
% extinit  external camera parameters for initialization
% cmodel   camera model with radial and tangential distortion
% frames   camera model for multiple views
% lmoptc   Levenberg-Marquardt optimization for cacal and circal
% jacobc   produces the Jacobian matrix for optimization
% cinit    forms the matrix containing the control point position, 
%          orientation an geometry
%
% Demonstration:
%
% caldemo  performs two example calibrations

% Author: Janne Heikkila, University of Oulu, Finland

