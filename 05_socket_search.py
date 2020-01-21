#!/usr/bin/python

import socket

def get_ip_status(ip,port):
	print "checking " + ip;
	server = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
	try:
		server.connect ((ip,port))
		print ('{0} port {1} is open'.format(ip,port))
	except Exception as err:
		print ('{0} port {1} is not open'.format(ip,port))
	finally:
		server.close ()

for i in range (0,256):
	ip = "172.16.43." + str(i)
	get_ip_status (ip, 22)
