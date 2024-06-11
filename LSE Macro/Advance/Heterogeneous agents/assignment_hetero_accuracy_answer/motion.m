%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************
function kp = motion(k,z)

%True law of motion
kp = max(0.2 + 0.8 * k + 0.18 * log(k) + z,0.7);

end

