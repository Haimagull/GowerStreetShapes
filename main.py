"""
Docstring for repo_gowerstreetshapes.GowerStreetShapes.main
"""
import subprocess
import sys

def welcome() :
    image = """
                                                                                                               
                             .@...@-     .%@@.  %@.                                                                 
                         @.     .       ..@. .=%.   ..%.                                                            
       . .            -%.    ..  .@%%.       :@ .. -+.  .%.         ...               .                       .     
        ...        @.     @...    @@..   .@  ....** .@.    %.%=    ..                .. .                 ...       
         .       @  %%. . .+=.     ..      .@ .      %  .%  .@...++@. .               . .                .          
                 . @..  .     .@..=@@.. *@     :@@.    .@..%   .*..%..+                                ..           
                 %%..@.   .%..              .@.   .@.     .@..@.. .=.@*.@                                           
                  .     @.           @%@@+     ..-@    @.   .%..%   .@ .@.@.                                        
                  .@.  %.       @ .        .%      @.%%..#     @.%    =. @  @.                                      
                    %..     *@.      %@      .%...%      .@.     %.@    @..@  %.                                    
                   ..    .+.    %..@                                #..   @  @ .%                                   
           .       @.   ..     ..% %               .               @..     @..% .@.        .                        
         ..        @    . %%.     %..           ....             .@...:.    % .. ..@.   ...               ...       
        .          %    @... ..   .           ##+%.              .  .@.@     : %    %@@..                 . ..      
      ..            =   . ..*.   %.         .:..%.                .. .      .  @   ..  ..@               . ...      
     .             @.   *  @%.    @         #.  %                         %.. # %   . @.%.%.                        
                    %..  #.       @.        @   @      =.@     ..%-.     .%.  @.@     @.:. .                        
                       %. @#      %%%      .@   .%..   %  .%     @..@           .     .  *..                        
                       .@. .# @@    .       . .@. %.   #+ %#@    ..@.@.              @. @ %.*.                      
                         .# @ .@@   .*       =.    @%   .@. .. %    .%. .%%%...  .@%. .%@.%  .@                     
                          %  -@.%.    @.      .@%    %     .@  .@%%.  .%@                 %. %.                     
                            @  @     %. -       ...   @.     %      @  .  ..%.  .%@ @@.   @ ..  %%                  
        .                   . %..=.  .@@.%     .  =..   @.   .%%.   ..@..       .%%@%%% ..   . %.%. @       ...     
        ..                  .             @. .. ... .@ . .@.     .@..    .@.         ..... ..   ...@  %.   ..       
       ...                ..               ...+...      :%. % .  ...@..-%.  .@%%%%@@@@*%=..     @.# .  % ..         
                        ..                    %.         .%  ..@.    .@#..               @.   %@.. %   %.           
                                              ..@%        .@%.   @           @%@@.     .@..    #%.    @..           
                                                  %.            @ .@.            .. ..@@@.   . .    %..             
                                                   %.  %.         %.  .  @%@..                  =@+.                
                                                        %           :%.  %=@. .% -%@+#@@%%.. .                      
                                                         .@@%@        .%@..#@%@@.  @%                               
                                                              ..  .%.  .%#%@@%%@                         
  ____                          ____  _                 _   
 / ___| _____      _____ _ __  / ___|| |_ _ __ ___  ___| |_ 
| |  _ / _ \ \ /\ / / _ \ '__| \___ \| __| '__/ _ \/ _ \ __|
| |_| | (_) \ V  V /  __/ |     ___) | |_| | |  __/  __/ |_ 
 \____|\___/ \_/\_/ \___|_|    |____/ \__|_|  \___|\___|\__|
/ ___|| |__   __ _ _ __   ___  ___                          
\___ \| '_ \ / _` | '_ \ / _ \/ __|                         
 ___) | | | | (_| | |_) |  __/\__ \                         
|____/|_| |_|\__,_| .__/ \___||___/                         
                  |_|                                       
    """
    print(image)
    print("Welcome to GowerStreetShapes. This tool is meant to help you navigate through fof outputs of Gower Street Simulations to study shapes. " \
    "If you need more info, take a look at https://github.com/Haimagull/GowerStreetShapes. Have fun ! - O.L.")

def decision() :
    options = """
    Here are the tools available :

    (1) haloshaper : visualize influence of semi-axis ratios parameters on an ellipsoid.
    (2) ratiosmoments : get main 1st order moments to visualise evolution across redshift, you can use it on shapez outputs.
    (3) shapedistribz : plot distribution of halo shape, you can use it on shapez outputs.
    (4) shapez : generates csv file with semi axis of halos at each redshift. Avoids great computing time when later trying to plot.
    (5) corrfunction : plot correlation functions from Covo output.
    (6) foftocovo : converts a FOF output into csv readable by Covo (3D correlation code, Kai Hoffman, on Bitbucket).
    (7) foftoIACorr : converts a FOF output into csv readable by IACorr (ellipticity correlation code using treecorr, elisabethjg, on Github).
    """
    print(options)
    number = int(input("Please, indicate what you want to use with its number: "))
    if number == 1 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/UnderstandingShapes/HaloShaper.py"])
    elif number == 2 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/HaloShapeDistribution/ratios_moments.py"])
    elif number == 3 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/HaloShapeDistribution/shape_distribution_z.py"])
    elif number == 4 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/HaloShapeDistribution/shape_z.py"])
    elif number == 5 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/ShapeCorrelations/correlation_function.py"])
    elif number == 6 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/ShapeCorrelations/fof_to_covo.py"])
    elif number == 7 :
        subprocess.run([sys.executable, "GowerStreetShapes/src/ShapeCorrelations/fof_to_IACorr.py"])

if __name__ == "__main__" :
    welcome()
    decision()



