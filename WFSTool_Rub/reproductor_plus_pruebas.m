player = reproductor_plus();

comMat = [1 1];
player.setProps('comMatrix', comMat);
player.setProps('audioFileName', 'C:\Users\Rubén\Music\Ariana Grande - Dangerous Woman\01 - Moonlight.mp3', 1);
player.setProps('getDelayFun', @() [0; 0], [1 1]);
player.setProps('getAttenFun', @() [1; 1], [1 1]);
player.setProps('getDelayFun', @() [0; 0], [1 2]);
player.setProps('getAttenFun', @() [1; 1], [1 2]);
player.setProps('frameDuration', 2);

player.setProps('device', 'Altavoces (Dispositivo de High Definition Audio)', 1);
player.setProps('device', 'Controlador primario de sonido', 2);

setup(player, [], []);
delays = player.testDelayBetweenDevices();

order.action = 'stop';
player.executeOrder(order);

order.action = 'play';
player.executeOrder(order);

