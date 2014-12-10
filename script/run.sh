#! /bin/bash

./bin/release/collision_handle ./data/dino/stair.ini  >  tempt/dino_stair.txt
./bin/release/collision_handle ./data/dino/plane.ini  >  tempt/dino_plane.txt
./bin/release/collision_handle ./data/dino/cylinder.ini  >  tempt/dino_cylinder.txt
./bin/release/collision_handle ./data/dino/ball.ini  >  tempt/dino_ball.txt

./bin/release/collision_handle ./data/dragon_asia/stair.ini  >  tempt/dragon_asia_stair.txt
./bin/release/collision_handle ./data/dragon_asia/plane.ini  >  tempt/dragon_asia_plane.txt

./bin/release/collision_handle ./data/beam/ball.ini  >  tempt/beam_ball.txt

./bin/release/collision_handle ./data/dragon/stair.ini  >  tempt/dragon_stair.txt
./bin/release/collision_handle ./data/dragon/bowl.ini  >  tempt/dragon_bowl.txt
./bin/release/collision_handle ./data/dragon/plane.ini  >  tempt/dragon_plane.txt
