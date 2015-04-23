import os

class ProdigalHit:
	def __init__(self,
				 start,
				 end,
				 strand,
				 point,
				 ):
		# Calculate distance to start
		if start <= p <= end:
			self.distance = 0
		else:
			self.distance = min(abs(start - p), abs(end - p))
		
				 
				  
