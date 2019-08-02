% this function indexes all ways to assign n balls to d levels such that each level has at least one ball and the first level always have one ball
% there are nchoosek( n-2 , d-2 ) many possible options for ix
% the ix variable can take values 1, 2, 3, 4, 5, ...

function [balls_levels] = index_n_balls_d_levels_one_ball_1st_atleastone_later(ix , numBalls , numLevels)
    
    n = numBalls  - 2;
    k = numLevels - 2;
    
    [comb] = index_n_choose_k(ix , n , k);  % comb is a vector
    
    balls_levels = [1, diff([0 , comb, numBalls-1])];  %balls_level is a vector
    
end