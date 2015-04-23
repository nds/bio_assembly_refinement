import os

class ProdigalHit:
	def __init__(self,
				 start,
				 end,
				 strand,
				 point,
				 ):
		self.start = start
		self.end = end
		self.strand = strand
		# Calculate distance to start
		if start <= point <= end:
			self.distance = 0
		else:
			self.distance = min(abs(start - point), abs(end - point))
		
				 
				  
