# RetroBioCat-DB
Please note - the code in this repository is **not** licensed for commercial use - academic use only is allowed.
This repository is an extension of our original RetroBioCat work - which still available at https://github.com/willfinnigan/RetroBioCat

This repository holds additional code pertaining to the database aspects of RetroBioCat, including data curation and analysis.

We recommend using retrobiocat through the online version hosted at https://retrobiocat.com  

## Installation

### 1. Copy necessary docker files  
Clone this repository and move to `/docker` directory.  

### 2. Pull the latest docker container for retrobiocat_db
Run `sudo docker pull willfinnigan/retrobiocat_db:latest`.

Alternatively, you can build your own version of this container using the provided Dockefile. 
You will need to edit the docker-compose.yml file to utilise your built container if you do this.  

### 3. Edit configuration as required
Edit docker_env.env as necessary to change the configuration.
`sudo nano docker_env.env`  
It is recommended to alter the `SECRET_PASSWORD` and `SECURITY_PASSWORD_SALT` from the default settings for security.

### 4. Run retrobiocat
Run `sudo docker-compose --env-file docker_env.env up`.  
RetroBioCat should now be available at the [localhost address](http://127.0.0.1) 
(if you have changed the NGINX_PORT variable this must also now be included)  

### 5. Initialise the database
- Login using the default admin account (setup in the config below).  
- At the right hand side of the nav bar, click the user and go to `Initialise Database` in the menu.  
- Load the `mongo_dump.gz` file. The file in this repo is only a trial version. Check the terminal output for successful database initialisation.
- You will now need to log in again.  
- Return to `Initialise Database`, click through the `Download additional required files` buttons to queue jobs to download these files. 
Again, check the terminal output for any errors.
- Restart all the docker containers (the downloaded files are not correctly loaded otherwise). RetroBioCat should now be ready for use.

## Config options
### Important options to change
Please edit the `SECRET_KEY` and `SECURITY_PASSWORD_SALT` variables.  
Note, changing the secret key and password salt once accounts have already been set up will stop existing passwords from working.  

### Admin account
The admin account which is created when retrobiocat is started.  Again, defaults will work but can be changed. 
Please note, once an account is created, changing settings here will simply create a new account without deleting the old one.  

`ADMIN_EMAIL=admin1@email.com` - admin email address (can be anything).  
`ADMIN_PASSWORD=password1` - admin password.





